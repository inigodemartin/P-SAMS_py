#!/usr/bin/env python3
"""
Tests for screening off-targets against more than one transcriptome.

The case they cover is the one Alberto described: a target gene belonging
to a pathogen that infects several hosts (a TSWV gene, say), where a guide
is only usable when it is clean in every host at once. Running the design
once per host and crossing the results by hand is what these replace.

Like test_pepper_bugs.py, everything is built from synthetic FASTA and
SQLite inputs in a temporary directory, so no transcriptome, annotation
database or network access is needed.

Run from the repository root with:  python3 test/test_multi_background.py
"""

import json
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.pipeline import (  # noqa: E402
    _load_tf_cache,
    _tf_cache_key,
    _tf_cache_store,
    _run_targetfinder,
    get_tsites,
    serial_jobs,
)
from src.utils import Background, parse_species  # noqa: E402

TARGETFINDER = REPO_ROOT / "TargetFinder" / "targetfinder.pl"

SEED = 15
SITE_LENGTH = 21
OFFSET = SITE_LENGTH - SEED - 1

# Same construction as test_pepper_bugs.py: a 21-nt guide and a transcript
# carrying its exact reverse complement, so TargetFinder is guaranteed a
# perfect hit without real sequence data.
GUIDE = "TTCGACGTACGTACGTACGTA"
_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
TARGET_SITE = "".join(_COMPLEMENT[base] for base in reversed(GUIDE))
TRANSCRIPT = "ACGTACGGTTACCGGTTACCG" * 3 + TARGET_SITE + "GGTTACCGGTTACCGGTTACC" * 3

# A transcript sharing no homology with the guide, so TargetFinder reports
# "No results" against it. That is what a host transcriptome looks like to
# a guide designed against a pathogen gene it does not carry.
UNRELATED = "TTTTCCCCAAAATTTTCCCCAAAATTTTCCCCAAAATTTTCCCCAAAATTTTCCCCAAAA"


def write_fasta(path: Path, entries) -> Path:
    """entries: {header: sequence}"""
    path.write_text("".join(f">{h}\n{seq}\n" for h, seq in entries.items()))
    return path


def make_annotation_db(path: Path, rows) -> sqlite3.Connection:
    """Build the minimal annotation table off_target_check() queries."""
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE annotation (transcript TEXT PRIMARY KEY, description TEXT)")
    conn.executemany("INSERT INTO annotation VALUES (?, ?)", rows)
    conn.commit()
    conn.row_factory = sqlite3.Row
    return conn


def make_kmer_db(path: Path, rows) -> sqlite3.Connection:
    """Build the minimal kmers table get_tsites() queries. rows: {kmer: transcripts}"""
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE kmers (kmer TEXT PRIMARY KEY, transcripts TEXT NOT NULL)")
    conn.executemany("INSERT INTO kmers VALUES (?, ?)", list(rows.items()))
    conn.commit()
    conn.row_factory = sqlite3.Row
    return conn


def seed_of(site: str) -> str:
    return site[OFFSET:OFFSET + SEED]


def run_candidate(tmp: Path, backgrounds, names="TX00001", ids=None):
    """Push one candidate through serial_jobs against the given backgrounds."""
    tf_dir = tmp / "tf_results"
    tf_dir.mkdir(exist_ok=True)
    (tmp / ".cache").mkdir(exist_ok=True)

    site = {
        "guide": GUIDE,
        "star": GUIDE,
        "oligo1": "oligo1",
        "oligo2": "oligo2",
        "names": names,
        "seqs": TARGET_SITE,
    }

    return serial_jobs(
        target_count=1,
        construct="amiRNA",
        ids=ids if ids is not None else {"TX00001": TRANSCRIPT},
        site_scores=[site],
        targetfinder=str(TARGETFINDER),
        backgrounds=backgrounds,
        bg=True,
        subopt_name=str(tmp / "subopt.tsv"),
        opt_name=str(tmp / "opt.tsv"),
        output_folder=str(tmp),
        accession_list=["TX00001"],
        tf_dir=str(tf_dir),
        potential_target_n=1,
        fasta=True,
    )


class SeedFilterTest(unittest.TestCase):
    """
    The k-mer filter is the first of the two off-target layers and runs
    before any guide is designed. With several hosts a seed has to be
    absent from all of them, not just from the first.
    """

    def _sites(self, tmp, backgrounds):
        return get_tsites({"viral_gene": TRANSCRIPT}, SEED, bg=True, backgrounds=backgrounds)

    def test_a_seed_present_only_in_the_second_host_is_still_discarded(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            # The seed of the very first window of the input sequence.
            first_site = TRANSCRIPT[0:SITE_LENGTH]
            burned = seed_of(first_site)

            empty = make_kmer_db(tmp / "a.db", {})
            carrier = make_kmer_db(tmp / "b.db", {burned: "HOST_B_TX1"})

            host_a = Background("host_a", empty, str(tmp / "a.fa"))
            host_b = Background("host_b", carrier, str(tmp / "b.fa"))

            against_first = self._sites(tmp, [host_a])
            against_both = self._sites(tmp, [host_a, host_b])

            empty.close()
            carrier.close()

            self.assertIn(first_site, [s["seq"] for s in against_first],
                          "the site survives when only the clean host is screened")
            self.assertNotIn(first_site, [s["seq"] for s in against_both],
                             "a seed found in the second host must be discarded too")

    def test_a_seed_belonging_to_a_target_is_kept(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            first_site = TRANSCRIPT[0:SITE_LENGTH]
            # The k-mer is registered, but against the input sequence
            # itself, so it is an on-target and must not be filtered out.
            own = make_kmer_db(tmp / "a.db", {seed_of(first_site): "viral_gene"})
            host = Background("host_a", own, str(tmp / "a.fa"))

            sites = self._sites(tmp, [host])
            own.close()

            self.assertIn(first_site, [s["seq"] for s in sites])


class OffTargetUnionTest(unittest.TestCase):
    """A candidate is optimal only when no host has an off-target for it."""

    def test_a_candidate_clean_in_one_host_but_not_the_other_is_suboptimal(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            db_a = write_fasta(tmp / "a.fa", {"TX00001": TRANSCRIPT})
            db_b = write_fasta(tmp / "b.fa", {"TXB0001": TRANSCRIPT})
            conn_a = make_annotation_db(tmp / "a.db", [("TX00001", "intended target")])
            conn_b = make_annotation_db(tmp / "b.db", [("TXB0001", "a pepper gene")])

            host_a = Background("tomato", conn_a, str(db_a))
            host_b = Background("pepper", conn_b, str(db_b))

            opt, subopt = run_candidate(tmp, [host_a, host_b])
            conn_a.close()
            conn_b.close()

            self.assertEqual(len(opt), 0, "the pepper off-target must disqualify it")
            self.assertEqual(len(subopt), 1)

            site = subopt[0]["site"]
            self.assertEqual(site["offtarget_number"], 1)
            self.assertEqual(site["offtarget_list"], ["pepper:TXB0001"],
                             "the off-target has to name the host it was found in")

            document = json.loads(site["tf"])
            hits = document[next(iter(document))]["hits"]
            species = {hit.get("Target species") for hit in hits if "Target species" in hit}
            self.assertIn("pepper", species)

    def test_the_same_candidate_is_optimal_against_the_clean_host_alone(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            db_a = write_fasta(tmp / "a.fa", {"TX00001": TRANSCRIPT})
            conn_a = make_annotation_db(tmp / "a.db", [("TX00001", "intended target")])
            host_a = Background("tomato", conn_a, str(db_a))

            opt, subopt = run_candidate(tmp, [host_a])
            conn_a.close()

            self.assertEqual(len(opt), 1, "with one host there is nothing to disqualify it")
            self.assertEqual(len(subopt), 0)

    def test_a_single_background_run_does_not_label_hits_with_a_species(self):
        """Screening against one transcriptome keeps the output it always had."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            db_a = write_fasta(tmp / "a.fa", {"TX00001": TRANSCRIPT, "TX00002": TRANSCRIPT})
            conn_a = make_annotation_db(
                tmp / "a.db", [("TX00001", "intended target"), ("TX00002", "off-target")]
            )
            host_a = Background("tomato", conn_a, str(db_a))

            opt, subopt = run_candidate(tmp, [host_a])
            conn_a.close()

            self.assertEqual(len(subopt), 1)
            site = subopt[0]["site"]
            self.assertEqual(site["offtarget_list"], ["TX00002"],
                             "a single-host run keeps bare accessions")

            document = json.loads(site["tf"])
            hits = document[next(iter(document))]["hits"]
            self.assertFalse([hit for hit in hits if "Target species" in hit])


class PathogenTargetTest(unittest.TestCase):
    """
    The target gene of a pathogen is in none of the host transcriptomes, so
    TargetFinder reporting nothing is the best possible outcome. Before
    multi-background support, that answer discarded the candidate instead,
    which silently threw away every perfectly clean guide.
    """

    def test_a_guide_absent_from_every_host_is_optimal(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            db_a = write_fasta(tmp / "a.fa", {"HOST_A_TX1": UNRELATED})
            db_b = write_fasta(tmp / "b.fa", {"HOST_B_TX1": UNRELATED})
            conn_a = make_annotation_db(tmp / "a.db", [("HOST_A_TX1", "a tomato gene")])
            conn_b = make_annotation_db(tmp / "b.db", [("HOST_B_TX1", "a pepper gene")])

            host_a = Background("tomato", conn_a, str(db_a))
            host_b = Background("pepper", conn_b, str(db_b))

            opt, subopt = run_candidate(
                tmp, [host_a, host_b],
                names="viral_gene", ids={"viral_gene": TRANSCRIPT},
            )
            conn_a.close()
            conn_b.close()

            self.assertEqual(len(subopt), 0)
            self.assertEqual(len(opt), 1, "a guide with no host off-target must be optimal")

            # Its on-target pairing is still reported, computed against the
            # input sequence rather than against any host transcriptome.
            document = json.loads(opt[0]["tf"])
            hits = document[next(iter(document))]["hits"]
            self.assertEqual(len(hits), 1)
            self.assertEqual(hits[0]["Target accession"], "viral_gene")

    def test_one_host_carrying_a_hit_still_disqualifies_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            db_a = write_fasta(tmp / "a.fa", {"HOST_A_TX1": UNRELATED})
            db_b = write_fasta(tmp / "b.fa", {"HOST_B_TX1": TRANSCRIPT})
            conn_a = make_annotation_db(tmp / "a.db", [("HOST_A_TX1", "a tomato gene")])
            conn_b = make_annotation_db(tmp / "b.db", [("HOST_B_TX1", "a pepper gene")])

            opt, subopt = run_candidate(
                tmp,
                [Background("tomato", conn_a, str(db_a)),
                 Background("pepper", conn_b, str(db_b))],
                names="viral_gene", ids={"viral_gene": TRANSCRIPT},
            )
            conn_a.close()
            conn_b.close()

            self.assertEqual(len(opt), 0)
            self.assertEqual(subopt[0]["site"]["offtarget_list"], ["pepper:HOST_B_TX1"])


class AccessionBranchTest(unittest.TestCase):
    """
    With -a the target is a transcript of the searched transcriptome, so
    that branch keeps its old meaning: a guide aligning nowhere does not
    even hit its own target and is dropped.
    """

    def _run(self, tmp, backgrounds, names="TX00001"):
        tf_dir = tmp / "tf_results"
        tf_dir.mkdir(exist_ok=True)
        (tmp / ".cache").mkdir(exist_ok=True)
        site = {
            "guide": GUIDE, "star": GUIDE, "oligo1": "o1", "oligo2": "o2",
            "names": names, "seqs": TARGET_SITE,
        }
        return serial_jobs(
            target_count=1, construct="amiRNA", ids={names: TRANSCRIPT},
            site_scores=[site], targetfinder=str(TARGETFINDER),
            backgrounds=backgrounds, bg=True,
            subopt_name=str(tmp / "subopt.tsv"), opt_name=str(tmp / "opt.tsv"),
            output_folder=str(tmp), accession_list=[names], tf_dir=str(tf_dir),
            potential_target_n=1, fasta=False,
        )

    def test_a_candidate_hitting_only_its_own_target_is_optimal(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            db = write_fasta(tmp / "a.fa", {"TX00001": TRANSCRIPT})
            conn = make_annotation_db(tmp / "a.db", [("TX00001", "intended target")])

            opt, subopt = self._run(tmp, [Background("tomato", conn, str(db))])
            conn.close()

            self.assertEqual(len(opt), 1)
            self.assertEqual(len(subopt), 0)
            document = json.loads(opt[0]["tf"])
            self.assertEqual(document[next(iter(document))]["hits"][0]["Target accession"], "TX00001")

    def test_a_guide_aligning_nowhere_is_still_dropped(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            db = write_fasta(tmp / "a.fa", {"TX00001": UNRELATED})
            conn = make_annotation_db(tmp / "a.db", [("TX00001", "intended target")])

            opt, subopt = self._run(tmp, [Background("tomato", conn, str(db))])
            conn.close()

            self.assertEqual(opt, [])
            self.assertEqual(subopt, [])


class CacheKeyTest(unittest.TestCase):
    """
    A TargetFinder result depends on the transcriptome it was searched
    against, so reusing it across hosts would report a candidate as clean
    in a host it was never checked in.
    """

    def test_a_cached_result_is_not_reused_for_another_background(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            cache = _load_tf_cache(tmp / "cache.json")
            _tf_cache_store(cache, "tomato", GUIDE, "amiRNA0", ["tomato result"])

            self.assertIn(_tf_cache_key("tomato", GUIDE), cache["results"])
            self.assertNotIn(_tf_cache_key("pepper", GUIDE), cache["results"])

    def test_a_cached_result_is_reused_for_its_own_background(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            cache = _load_tf_cache(tmp / "cache.json")
            _tf_cache_store(cache, "tomato", GUIDE, "placeholder", ['{"placeholder": {}}'])

            # A missing targetfinder path would raise if it were invoked,
            # so returning the relabelled cached lines proves it was not.
            lines = _run_targetfinder("/nonexistent/targetfinder.pl", GUIDE,
                                      Background("tomato", None, "/nonexistent.fa"),
                                      "amiRNA0", cache)
            self.assertEqual(lines, ['{"amiRNA0": {}}'])

    def test_the_cache_survives_a_round_trip_through_disk(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            cache = _load_tf_cache(tmp / "cache.json")
            _tf_cache_store(cache, "tomato", GUIDE, "amiRNA0", ["tomato result"])

            reloaded = _load_tf_cache(tmp / "cache.json")
            self.assertEqual(reloaded["results"][_tf_cache_key("tomato", GUIDE)],
                             ("amiRNA0", ["tomato result"]))

    def test_a_cache_written_before_multi_background_support_is_dropped(self):
        """
        Old entries are keyed by guide alone and carry no record of the
        transcriptome they were computed against, so they cannot be trusted.
        """
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            path = tmp / "cache.json"
            path.write_text(json.dumps({GUIDE: {"query_name": "amiRNA0", "lines": ["old"]}}))

            self.assertEqual(_load_tf_cache(path)["results"], {})


class ConfigValidationTest(unittest.TestCase):
    """
    Every species screened against needs its own database, built by
    create_db.py and registered in psams.conf. sqlite3 creates an empty
    file for a path that does not exist, so a wrong path used to surface
    only once the run was well under way, as a bare "no such table".
    """

    def _config(self, body):
        import configparser
        config = configparser.ConfigParser()
        config.read_string(body)
        return config

    def _complete_species(self, tmp, name):
        fasta = write_fasta(tmp / f"{name}.fa", {"TX00001": TRANSCRIPT})
        db = tmp / f"{name}.db"
        conn = sqlite3.connect(db)
        conn.execute("CREATE TABLE kmers (kmer TEXT PRIMARY KEY, transcripts TEXT NOT NULL)")
        conn.execute("CREATE TABLE annotation (transcript TEXT PRIMARY KEY, description TEXT)")
        conn.commit()
        conn.close()
        return f"[{name}]\nmRNA={fasta}\nsql={db}\n"

    def test_a_registered_species_connects(self):
        from src.utils import connect_backgrounds
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            config = self._config(self._complete_species(tmp, "tomato"))
            backgrounds = connect_backgrounds(config, "tomato")

            self.assertEqual([b.name for b in backgrounds], ["tomato"])
            backgrounds[0].conn.close()

    def test_an_unregistered_species_is_refused(self):
        from src.utils import connect_backgrounds
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            config = self._config(self._complete_species(tmp, "tomato"))
            with self.assertRaises(SystemExit):
                connect_backgrounds(config, "tomato,pepper")

    def test_a_missing_database_is_refused_before_the_run_starts(self):
        from src.utils import connect_backgrounds
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            fasta = write_fasta(tmp / "pepper.fa", {"TX00001": TRANSCRIPT})
            config = self._config(
                self._complete_species(tmp, "tomato")
                + f"\n[pepper]\nmRNA={fasta}\nsql={tmp / 'absent.db'}\n"
            )
            with self.assertRaises(SystemExit):
                connect_backgrounds(config, "tomato,pepper")

            self.assertFalse((tmp / "absent.db").exists(),
                             "the missing path must not be created as an empty database")

    def test_a_database_without_the_expected_tables_is_refused(self):
        from src.utils import connect_backgrounds
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            fasta = write_fasta(tmp / "pepper.fa", {"TX00001": TRANSCRIPT})
            empty = tmp / "empty.db"
            sqlite3.connect(empty).close()
            config = self._config(f"[pepper]\nmRNA={fasta}\nsql={empty}\n")

            with self.assertRaises(SystemExit):
                connect_backgrounds(config, "pepper")

    def test_a_missing_transcript_fasta_is_refused(self):
        from src.utils import connect_backgrounds
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            body = self._complete_species(tmp, "tomato")
            body = body.replace(str(tmp / "tomato.fa"), str(tmp / "absent.fa"))
            with self.assertRaises(SystemExit):
                connect_backgrounds(self._config(body), "tomato")


class ArgumentTest(unittest.TestCase):
    """-s/--species takes a list, but only where the target is unambiguous."""

    def _parse(self, argv):
        from src.args_parser import parse_args
        saved = sys.argv
        sys.argv = ["psams.py"] + argv
        try:
            return parse_args()
        finally:
            sys.argv = saved

    def test_species_list_is_split(self):
        self.assertEqual(parse_species("tomato, pepper"), ["tomato", "pepper"])
        self.assertEqual(parse_species(""), [])
        self.assertEqual(parse_species(None), [])

    def test_fasta_input_accepts_several_species(self):
        args = self._parse(["-f", "gene.fa", "-s", "tomato,pepper"])
        self.assertEqual(parse_species(args.species), ["tomato", "pepper"])

    def test_accessions_refuse_several_species(self):
        with self.assertRaises(SystemExit):
            self._parse(["-a", "GENE1", "-s", "tomato,pepper"])

    def test_disabling_offtargets_refuses_several_species(self):
        with self.assertRaises(SystemExit):
            self._parse(["-f", "gene.fa", "-s", "tomato,pepper", "-n"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
