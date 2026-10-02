#!/usr/bin/env python3
"""
Regression tests for the three bugs found while running P-SAMS against the
pepper (Capsicum annuum) transcriptome.

Each test fails against the code as it was before the fix:

  1. An unquoted database path reaching the shell. A path containing spaces
     or parentheses (a Google Drive folder, for instance) broke the
     ssearch36 call silently and the pipeline reported a false
     "0 optimal sites".
  2. Unescaped fields in TargetFinder's JSON output. FASTA headers holding a
     backslash (two pepper genes do) produced invalid JSON. The same bug has
     a second, independent home in off_target_check(), on the Python side.
  3. The -f plus suboptimal-candidate branch of serial_jobs() storing the
     TargetFinder result as a list instead of joined text, which raised a
     TypeError whenever off-targets were checked.
  4. A crash leaving a 0-byte tf_results/*.json behind, which the resume path
     then loaded as an empty string and carried on with, making the gap
     permanent across reruns.

Run from the repository root with:  python3 test/test_pepper_bugs.py
"""

import json
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

TARGETFINDER = REPO_ROOT / "TargetFinder" / "targetfinder.pl"

# A 21-nt guide and a transcript carrying its exact reverse complement, so
# TargetFinder is guaranteed a perfect hit without needing real sequence data.
GUIDE = "TTCGACGTACGTACGTACGTA"
_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
TARGET_SITE = "".join(_COMPLEMENT[base] for base in reversed(GUIDE))
TRANSCRIPT = "ACGTACGGTTACCGGTTACCG" * 3 + TARGET_SITE + "GGTTACCGGTTACCGGTTACC" * 3


def write_db(path: Path, headers) -> Path:
    """Write a small FASTA database, one transcript per supplied header."""
    path.write_text("".join(f">{h}\n{TRANSCRIPT}\n" for h in headers))
    return path


def run_targetfinder(db: Path) -> str:
    """Run targetfinder.pl in JSON mode against db and return its stdout."""
    result = subprocess.run(
        [
            "perl", str(TARGETFINDER),
            "-s", GUIDE,
            "-d", str(db),
            "-q", "testquery",
            "-p", "json",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.stdout


def make_annotation_db(path: Path, rows) -> Path:
    """Build the minimal annotation table off_target_check() queries."""
    conn = sqlite3.connect(path)
    conn.execute(
        "CREATE TABLE annotation (transcript TEXT PRIMARY KEY, description TEXT)"
    )
    conn.executemany("INSERT INTO annotation VALUES (?, ?)", rows)
    conn.commit()
    conn.close()
    return path


class ShellQuotingTest(unittest.TestCase):
    """Bug 1 — database paths with spaces or parentheses."""

    def test_database_path_with_spaces_and_parentheses(self):
        with tempfile.TemporaryDirectory() as tmp:
            # Mimics the shape of a synced Google Drive folder.
            awkward = Path(tmp) / "My Drive (copy)"
            awkward.mkdir()
            db = write_db(awkward / "db.fa", ["TX00001 a normal gene"])

            stdout = run_targetfinder(db)

            self.assertNotIn(
                "No results for", stdout,
                "the database path was passed to the shell unquoted, so "
                "ssearch36 never ran and the hit was silently lost",
            )
            self.assertIn("TX00001", json.loads(stdout)["testquery"]["hits"][0]["Target accession"])


class JsonEscapingTest(unittest.TestCase):
    """Bug 2 — backslashes and quotes in FASTA headers and descriptions."""

    def test_targetfinder_json_survives_backslash_in_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            db = write_db(
                Path(tmp) / "db.fa",
                ['TX00001 gene with a back\\slash and a "quote"'],
            )

            stdout = run_targetfinder(db)

            # The point of the test: the output has to be parseable at all.
            parsed = json.loads(stdout)
            self.assertEqual(
                parsed["testquery"]["hits"][0]["Target accession"],
                'TX00001 gene with a back\\slash and a "quote"',
            )

    def test_off_target_check_json_survives_backslash_in_description(self):
        from src.utils import off_target_check

        with tempfile.TemporaryDirectory() as tmp:
            db = make_annotation_db(
                Path(tmp) / "annot.db",
                [("TX00001", 'kinase, back\\slash and a "quote"')],
            )
            conn = sqlite3.connect(db)
            conn.row_factory = sqlite3.Row

            tf_results = json.dumps(
                {
                    "amiRNA1": {
                        "hits": [
                            {
                                "Target accession": "TX00001",
                                "Score": "0",
                                "Coordinates": "64-84",
                                "Strand": "+",
                                "Target sequence": TARGET_SITE,
                                "Base pairing": ":" * 21,
                                "amiRNA sequence": GUIDE,
                            }
                        ]
                    }
                },
                indent=4,
            ).splitlines(keepends=True)

            site = {"names": "TX00001", "seqs": TARGET_SITE, "guide": GUIDE}
            _off, _on, json_lines, _list = off_target_check(
                site, tf_results, conn, "amiRNA"
            )
            conn.close()

            parsed = json.loads("\n".join(json_lines))
            hit = parsed["amiRNA1"]["hits"][0]
            self.assertEqual(hit["Target description"], 'kinase, back\\slash and a "quote"')


class FastaSuboptimalBranchTest(unittest.TestCase):
    """Bug 3 — -f together with off-target checking on a suboptimal candidate."""

    def test_fasta_mode_suboptimal_candidate_writes_valid_json(self):
        from src.pipeline import serial_jobs

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            # Two identical transcripts: one is the intended target, the
            # other is an unavoidable off-target, which forces the candidate
            # down the suboptimal branch.
            db = write_db(tmp / "db.fa", ["TX00001", "TX00002"])
            annot = make_annotation_db(
                tmp / "annot.db",
                [("TX00001", "intended target"), ("TX00002", "off-target")],
            )
            conn = sqlite3.connect(annot)
            conn.row_factory = sqlite3.Row

            tf_dir = tmp / "tf_results"
            tf_dir.mkdir()

            site = {
                "guide": GUIDE,
                "star": GUIDE,
                "oligo1": "oligo1",
                "oligo2": "oligo2",
                "names": "TX00001",
                "seqs": TARGET_SITE,
            }

            opt, subopt = serial_jobs(
                target_count=1,
                construct="amiRNA",
                ids={"TX00001": TRANSCRIPT},
                site_scores=[site],
                targetfinder=str(TARGETFINDER),
                mRNA_fa=str(db),
                conn=conn,
                bg=True,
                subopt_name=str(tmp / "subopt.tsv"),
                opt_name=str(tmp / "opt.tsv"),
                output_folder=str(tmp),
                accession_list=["TX00001"],
                tf_dir=str(tf_dir),
                potential_target_n=1,
                fasta=True,
            )
            conn.close()

            self.assertEqual(len(subopt), 1, "the candidate should be suboptimal")
            # Before the fix this was a list, and json.loads() raised a
            # TypeError inside serial_jobs() before reaching this point.
            self.assertIsInstance(subopt[0]["site"]["tf"], str)
            json.loads(subopt[0]["site"]["tf"])

            cached = sorted(tf_dir.glob("*_TargetFinder_result.json"))
            self.assertTrue(cached, "no TargetFinder result file was written")
            json.loads(cached[0].read_text())


class ResumeAfterCrashTest(unittest.TestCase):
    """
    A candidate whose cached TargetFinder result is unreadable has to be
    evaluated again, and its stale row has to leave the checkpoint TSV: the
    TSVs are reopened in append mode, so a row left behind would list the
    same guide twice, once pointing at a cache that cannot be read.
    """

    OPT_HEADER = "Site_index\tGuide\tStar\tOligo1\tOligo2\tNames\tSeqs\n"
    SUB_HEADER = ("Site_index\tOfftarget_N\tOfftarget_list\tGuide\tStar\t"
                  "Oligo1\tOligo2\tNames\tSeqs\n")

    def _run(self, tmp, cached):
        """cached: {site index: file contents, or None to write no file}."""
        tf_dir = Path(tmp) / "tf_results"
        tf_dir.mkdir()
        opt = Path(tmp) / "opt.tsv"
        sub = Path(tmp) / "sub.tsv"

        rows = "".join(
            f"{i}\tGUIDE{i:04d}\tSTAR{i:04d}\to1\to2\tn\ts\n" for i in cached)
        opt.write_text(self.OPT_HEADER + rows)
        sub.write_text(self.SUB_HEADER)

        for i, text in cached.items():
            if text is not None:
                (tf_dir / f"site_{i:04d}_TargetFinder_result.json").write_text(text)

        from src.utils import load_resume_state
        return load_resume_state(str(tf_dir), str(opt), str(sub)), opt

    def test_an_empty_cache_sends_its_candidate_back_to_targetfinder(self):
        with tempfile.TemporaryDirectory() as tmp:
            (opt_sites, _, seen, start), _ = self._run(
                tmp, {1: '{"ok": 1}', 2: "", 3: '{"ok": 3}'})

            self.assertEqual([s["guide"] for s in opt_sites],
                             ["GUIDE0001", "GUIDE0003"])
            # Not marked as seen, so serial_jobs() leaves it in site_scores.
            self.assertNotIn("GUIDE0002", seen)
            # Numbering still continues past every row the previous run wrote.
            self.assertEqual(start, 3)

    def test_a_truncated_cache_counts_as_unusable_too(self):
        with tempfile.TemporaryDirectory() as tmp:
            (opt_sites, _, seen, _), _ = self._run(
                tmp, {1: '{"hits": [', 2: '{"ok": 2}'})
            self.assertEqual([s["guide"] for s in opt_sites], ["GUIDE0002"])
            self.assertNotIn("GUIDE0001", seen)

    def test_the_stale_row_is_removed_from_the_checkpoint_tsv(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, opt = self._run(tmp, {1: '{"ok": 1}', 2: "", 3: '{"ok": 3}'})

            lines = opt.read_text().splitlines()
            self.assertEqual(len(lines), 3, opt.read_text())
            self.assertNotIn("GUIDE0002", opt.read_text())
            # The rows that were fine are untouched, header included.
            self.assertTrue(lines[0].startswith("Site_index"))
            self.assertIn("GUIDE0001", lines[1])
            self.assertIn("GUIDE0003", lines[2])

    def test_a_clean_resume_rewrites_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, opt = self._run(tmp, {1: '{"ok": 1}', 2: '{"ok": 2}'})
            self.assertEqual(len(opt.read_text().splitlines()), 3)

    def test_a_missing_cache_file_is_treated_the_same_way(self):
        with tempfile.TemporaryDirectory() as tmp:
            (opt_sites, _, seen, _), opt = self._run(
                tmp, {1: '{"ok": 1}', 2: None})
            self.assertEqual([s["guide"] for s in opt_sites], ["GUIDE0001"])
            self.assertNotIn("GUIDE0002", seen)
            self.assertNotIn("GUIDE0002", opt.read_text())


if __name__ == "__main__":
    unittest.main(verbosity=2)
