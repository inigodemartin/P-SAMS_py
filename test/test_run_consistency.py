#!/usr/bin/env python3
"""
Regression tests for a run's parameters staying with the results they
produced.

An output folder is named after the gene(s) or -r/--run_name, and the run
key adds only the construct, so everything else about the command line is
forgotten once a run finishes. Two parameters change what the results mean:

  -t/--foldback picks the precursor an amiRNA is built into, eudicot
  (AtMIR390a) or monocot (OsMIR390). The guide and its TargetFinder result
  are the same either way; the amiRNA*, the hairpin and the vector the pair
  can be cloned into are not.

  -s/--species picks the background transcriptome(s) screened. A guide is
  only optimal when it has no off-target in them, so results from other
  hosts are a different search.

Neither was checked anywhere, so the same choice ended up stated twice and
could disagree:

  1. The -V vector menu. A vector carries one precursor, but any vector
     could be picked for any foldback: `-t eudicot` (the default) plus
     pMDC32B-OsMIR390-B/c produced the monocot overhangs around a eudicot
     star, in a file named after the monocot vector.
  2. The checkpoint's Star column. Resuming with the other -t rebuilt the
     oligos from the guide (as of 2feb357) but read the star back from the
     checkpoint, so one row carried the old foldback's amiRNA* next to the
     new foldback's oligos.
  3. The cached run. A second run of the same gene with the other foldback,
     or against another host, was answered with the first run's results —
     "already executed", or its cloning oligos rebuilt for a vector its
     star was never designed for.

Run from the repository root with:  python3 test/test_run_consistency.py
"""

import json
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.oligo_design import (  # noqa: E402
    AMIRNA_VECTORS, check_vector_foldback, make_amirna_oligos, select_vector)
from src.utils import (  # noqa: E402
    amirna_json, apply_vector_to_amirna_output, check_cached_run,
    create_opt_subopt, create_outputs, load_resume_state, oligo_designer,
    write_amirna_tsv)
from psams import _warn_stale_vector_files  # noqa: E402

# Same 21-nt guide the other suites use.
GUIDE = "TTCGACGTACGTACGTACGTA"

EUDICOT_STAR, EUDICOT_FWD, EUDICOT_REV = oligo_designer(GUIDE, "eudicot")
MONOCOT_STAR, MONOCOT_FWD, MONOCOT_REV = oligo_designer(GUIDE, "monocot")


class VectorAndFoldbackAgreeTest(unittest.TestCase):
    """Bug 1 — the vector and -t/--foldback are the same choice, twice."""

    def test_every_amirna_vector_declares_the_foldback_it_carries(self):
        for vector, spec in AMIRNA_VECTORS.items():
            self.assertIn(spec[3], ("eudicot", "monocot"), vector)

    def test_the_two_natural_loop_vectors_declare_the_right_one(self):
        """
        pMDC32B-AtMIR390a-B/c and pMDC32B-OsMIR390-B/c reuse their
        precursor's loop verbatim, so for each of them exactly one
        foldback reproduces oligo_designer()'s own pair — which is the
        foldback the table has to name.
        """
        for vector, expected in (("pMDC32B-AtMIR390a-B/c", "eudicot"),
                                 ("pMDC32B-OsMIR390-B/c", "monocot")):
            matching = []
            for foldback in ("eudicot", "monocot"):
                star, fwd, rev = oligo_designer(GUIDE, foldback)
                if make_amirna_oligos(GUIDE, star, vector) == (fwd, rev):
                    matching.append(foldback)
            self.assertEqual(matching, [expected], vector)
            self.assertEqual(AMIRNA_VECTORS[vector][3], expected)

    def test_the_chimera_a_mismatched_pair_used_to_produce(self):
        """
        The concrete discrepancy: the eudicot star in the monocot vector
        differs from the real monocot pair by the star's last paired base
        (A instead of C), which is exactly the position the OsMIR390
        foldback fixes.
        """
        chimera, _ = make_amirna_oligos(
            GUIDE, EUDICOT_STAR, "pMDC32B-OsMIR390-B/c")
        self.assertEqual(
            chimera,
            "CTTGTTCGACGTACGTACGTACGTAATGATGATCACATTCGTTATCTATTTTTTTACGTACGTAAGTACGTCGAA")
        self.assertEqual(
            MONOCOT_FWD,
            "CTTGTTCGACGTACGTACGTACGTAATGATGATCACATTCGTTATCTATTTTTTTACGTACGTAAGTACGTCGAC")
        differing = [i for i, (a, b) in enumerate(zip(chimera, MONOCOT_FWD)) if a != b]
        self.assertEqual(differing, [74])

    def test_a_vector_of_the_other_foldback_is_refused(self):
        with self.assertRaises(SystemExit) as caught:
            check_vector_foldback("pMDC32B-OsMIR390-B/c", "eudicot")
        message = str(caught.exception)
        self.assertIn("monocot", message)
        self.assertIn("pMDC32B-AtMIR390a-B/c", message)

        for vector in ("pMDC32B-AtMIR390a-B/c", "pMDC32B-BS-AtMIR390a-B/c",
                       "pMDC32B-BS-AtMIR390a-A18G-B/c"):
            with self.assertRaises(SystemExit, msg=vector):
                check_vector_foldback(vector, "monocot")

    def test_a_matching_vector_is_accepted(self):
        for vector, spec in AMIRNA_VECTORS.items():
            self.assertIsNone(check_vector_foldback(vector, spec[3]), vector)

    def test_a_cached_run_is_not_given_a_vector_of_the_other_foldback(self):
        """
        Same check on the path clone_vector.py and a -V rerun take, where
        the foldback comes from the cached run rather than from -t.
        """
        data = {
            "foldback": "eudicot",
            "optimal": {"amiRNA Optimal Result 1": {
                "amiRNA": GUIDE, "amiRNA*": EUDICOT_STAR,
                "Forward Oligo": EUDICOT_FWD, "Reverse Oligo": EUDICOT_REV,
                "TargetFinder": []}},
        }
        with self.assertRaises(SystemExit):
            apply_vector_to_amirna_output(data, "pMDC32B-OsMIR390-B/c")
        # Refused before anything was rewritten.
        self.assertEqual(
            data["optimal"]["amiRNA Optimal Result 1"]["Forward Oligo"], EUDICOT_FWD)
        self.assertNotIn("vector", data)

        apply_vector_to_amirna_output(data, "pMDC32B-BS-AtMIR390a-B/c")
        self.assertEqual(data["vector"], "pMDC32B-BS-AtMIR390a-B/c")
        self.assertEqual(
            len(data["optimal"]["amiRNA Optimal Result 1"]["Forward Oligo"]), 58)

    def test_a_cache_from_before_the_foldback_was_recorded_is_eudicot(self):
        """Same convention as insert_design.check_foldback."""
        def cached():
            return {"optimal": {"amiRNA Optimal Result 1": {
                "amiRNA": GUIDE, "amiRNA*": EUDICOT_STAR,
                "Forward Oligo": "", "Reverse Oligo": "", "TargetFinder": []}}}

        apply_vector_to_amirna_output(cached(), "pMDC32B-AtMIR390a-B/c")
        with self.assertRaises(SystemExit):
            apply_vector_to_amirna_output(cached(), "pMDC32B-OsMIR390-B/c")


class ResumeRebuildsTheStarTest(unittest.TestCase):
    """Bug 2 — the checkpoint's Star column outliving the foldback."""

    def _checkpoint(self, tmp, star):
        tf_dir = Path(tmp) / "tf_results"
        tf_dir.mkdir()
        (tf_dir / "site_0001_TargetFinder_result.json").write_text('{"ok": 1}')
        (tf_dir / "site_0002_TargetFinder_result.json").write_text('{"ok": 2}')
        subopt_name, opt_name = create_outputs("job_amiRNA", tmp)
        with open(opt_name, "a") as out:
            out.write(f"1\t{GUIDE}\t{star}\tTX00001\tSEQ\n")
        with open(subopt_name, "a") as out:
            out.write(f"2\t1\tTX00002\t{GUIDE}\t{star}\tTX00001\tSEQ\n")
        return tf_dir, opt_name, subopt_name

    def test_resuming_under_the_other_foldback_rebuilds_the_star(self):
        with tempfile.TemporaryDirectory() as tmp:
            tf_dir, opt_name, subopt_name = self._checkpoint(tmp, EUDICOT_STAR)
            opt, subopt, _, _ = load_resume_state(
                str(tf_dir), opt_name, subopt_name, foldback="monocot")

        for site in (opt[0], subopt[0]["site"]):
            # Before the fix this was EUDICOT_STAR, next to monocot oligos.
            self.assertEqual(site["star"], MONOCOT_STAR)
            self.assertEqual(site["oligo1"], MONOCOT_FWD)
            # The star reported and the star inside the oligo are the same
            # molecule: the pair carries it minus its "CA" 3' overhang.
            self.assertTrue(site["oligo1"].endswith(site["star"][:19]))

    def test_the_resumed_result_matches_a_fresh_run_of_that_foldback(self):
        """
        End to end on the artifact the user reads: amiRNA* and the oligo
        pair of one result must describe one construct.
        """
        with tempfile.TemporaryDirectory() as tmp:
            tf_dir, opt_name, subopt_name = self._checkpoint(tmp, EUDICOT_STAR)
            opt, subopt, _, _ = load_resume_state(
                str(tf_dir), opt_name, subopt_name, foldback="monocot")
            counts = create_opt_subopt(subopt, opt, "amiRNA")
            # A vector is named because that is now the only kind of
            # amiRNA result that carries cloning oligos at all (see
            # amirna_json); the monocot vector is the one this foldback
            # can be cloned into.
            data = amirna_json(counts[0], counts[1], counts[2], counts[3],
                               Path(tmp) / "out.json", foldback="monocot",
                               vector="pMDC32B-OsMIR390-B/c")

        entry = data["optimal"]["amiRNA Optimal Result 1"]
        self.assertEqual(entry["amiRNA*"], MONOCOT_STAR)
        self.assertEqual(entry["Forward Oligo"], MONOCOT_FWD)
        self.assertEqual(entry["Reverse Oligo"], MONOCOT_REV)

    def test_resuming_under_the_same_foldback_is_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            tf_dir, opt_name, subopt_name = self._checkpoint(tmp, EUDICOT_STAR)
            opt, _, _, _ = load_resume_state(
                str(tf_dir), opt_name, subopt_name, foldback="eudicot")
        self.assertEqual(opt[0]["star"], EUDICOT_STAR)
        self.assertEqual(opt[0]["oligo1"], EUDICOT_FWD)

    def test_syntasirna_keeps_the_star_it_stored(self):
        """No foldback is passed for syntasiRNA, which has no hairpin."""
        with tempfile.TemporaryDirectory() as tmp:
            tf_dir, opt_name, subopt_name = self._checkpoint(tmp, "STARSEQ")
            opt, subopt, _, _ = load_resume_state(
                str(tf_dir), opt_name, subopt_name)
        self.assertEqual(opt[0]["star"], "STARSEQ")
        self.assertEqual(opt[0]["oligo1"], "")
        self.assertEqual(subopt[0]["site"]["star"], "STARSEQ")


class CachedRunIsTiedToItsParametersTest(unittest.TestCase):
    """
    Bug 3 — one folder, one cache, two different questions.

    The folder is named after the gene (or -r/--run_name) and the run key
    adds only the construct, so a second run of the same input reuses or
    extends the first one's results whatever else changed on the command
    line. Two things change the answer: the foldback, and the background
    transcriptome(s) the candidates were screened against.
    """

    def _cache(self, tmp, foldback="eudicot", species=None, name="run_amiRNA_psams.json"):
        path = Path(tmp) / name
        star = EUDICOT_STAR if foldback == "eudicot" else MONOCOT_STAR
        amirna_json(
            1, 0,
            {1: {"guide": GUIDE, "star": star, "oligo1": "", "oligo2": "", "tf": "{}"}},
            {}, path, foldback=foldback, species=species)
        return path

    def test_a_eudicot_run_does_not_answer_a_monocot_request(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._cache(tmp, foldback="eudicot")
            with self.assertRaises(SystemExit) as caught:
                check_cached_run(path, foldback="monocot")
        message = str(caught.exception)
        self.assertIn("eudicot", message)
        self.assertIn("monocot", message)
        self.assertIn("-r/--run_name", message)

    def test_a_run_screened_against_another_host_is_not_reused(self):
        """
        The reproduced case: a folder whose guides were screened against
        'alpha' answering a run that asked for 'beta', whose optimal sites
        are not the same ones — a guide clean in alpha can have an
        off-target in beta.
        """
        with tempfile.TemporaryDirectory() as tmp:
            path = self._cache(tmp, species=["alpha"])
            with self.assertRaises(SystemExit) as caught:
                check_cached_run(path, foldback="eudicot", species=["beta"])
        message = str(caught.exception)
        self.assertIn("alpha", message)
        self.assertIn("beta", message)

    def test_dropping_or_adding_a_host_is_not_reused_either(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._cache(tmp, species=["alpha", "beta"])
            with self.assertRaises(SystemExit):
                check_cached_run(path, foldback="eudicot", species=["alpha"])
            # Same hosts in another order is the same search.
            self.assertIsNone(
                check_cached_run(path, foldback="eudicot", species=["beta", "alpha"]))

    def test_the_same_parameters_go_through(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._cache(tmp, foldback="monocot", species=["alpha"])
            self.assertIsNone(
                check_cached_run(path, foldback="monocot", species=["alpha"]))

    def test_syntasirna_is_not_checked_against_the_foldback(self):
        """Its results are guides only, the same under either foldback."""
        with tempfile.TemporaryDirectory() as tmp:
            path = self._cache(tmp, foldback="eudicot", species=["alpha"])
            self.assertIsNone(check_cached_run(path, foldback=None, species=["alpha"]))

    def test_nothing_cached_yet_goes_through(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertIsNone(
                check_cached_run(Path(tmp) / "absent.json", foldback="monocot"))
            self.assertIsNone(check_cached_run(None, foldback="monocot"))

    def test_a_cache_without_the_fields_is_let_through_as_before(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "run_amiRNA_psams.json"
            path.write_text(json.dumps({"optimal": {}}))
            # No "foldback" key: taken at the default, eudicot.
            self.assertIsNone(check_cached_run(path, foldback="eudicot"))
            with self.assertRaises(SystemExit):
                check_cached_run(path, foldback="monocot")
            # No "species" key: nothing to compare, so it still resumes.
            self.assertIsNone(
                check_cached_run(path, foldback="eudicot", species=["alpha"]))

    def test_the_run_json_records_what_it_was_screened_against(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self._cache(tmp, species=["alpha", "beta"])
            data = json.loads(path.read_text())
        self.assertEqual(data["species"], ["alpha", "beta"])
        self.assertEqual(data["foldback"], "eudicot")


def _one_result_json(tmp, name, **kwargs):
    """Write a one-result amiRNA JSON through the real builder."""
    site = {1: {"guide": GUIDE, "star": EUDICOT_STAR,
                "oligo1": EUDICOT_FWD, "oligo2": EUDICOT_REV, "tf": []}}
    return amirna_json(1, 0, site, {}, Path(tmp) / name,
                       no_offtarget=True, foldback="eudicot", **kwargs)


class OligosOnlyWhereTheVectorIsNamedTest(unittest.TestCase):
    """
    A result without --vector used to carry Forward/Reverse Oligo anyway,
    holding the bare foldback hairpin. That pair is not vector-neutral: it
    is exactly one vector's pair, under a filename and column headings that
    name none, beside the vector file's own pair of a different length.
    Same trap as the checkpoint's (2feb357), one file over.
    """

    def test_the_pair_left_out_was_one_specific_vector_all_along(self):
        """Why the vector-less pair could not stay: it was never neutral."""
        for foldback, vector in (("eudicot", "pMDC32B-AtMIR390a-B/c"),
                                 ("monocot", "pMDC32B-OsMIR390-B/c")):
            star, fwd, rev = oligo_designer(GUIDE, foldback)
            self.assertEqual((fwd, rev), make_amirna_oligos(GUIDE, star, vector),
                             f"{foldback}: the bare hairpin is {vector}'s pair")

    def test_no_vector_no_oligos_in_the_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            data = _one_result_json(tmp, "out.json")
        entry = data["results"]["amiRNA Result 1"]
        self.assertEqual(entry["amiRNA"], GUIDE)
        self.assertEqual(entry["amiRNA*"], EUDICOT_STAR)
        self.assertNotIn("Forward Oligo", entry)
        self.assertNotIn("Reverse Oligo", entry)

    def test_a_named_vector_still_carries_its_pair(self):
        with tempfile.TemporaryDirectory() as tmp:
            data = _one_result_json(tmp, "out.json",
                                    vector="pMDC32B-AtMIR390a-B/c")
        entry = data["results"]["amiRNA Result 1"]
        self.assertEqual(entry["Forward Oligo"], EUDICOT_FWD)
        self.assertEqual(entry["Reverse Oligo"], EUDICOT_REV)

    def test_the_tsv_drops_the_columns_with_the_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            tsv = Path(tmp) / "out.tsv"
            write_amirna_tsv(_one_result_json(tmp, "out.json"), tsv)
            header = tsv.read_text().splitlines()[0].split("\t")
            rows = tsv.read_text().splitlines()[1:]
        self.assertNotIn("Forward Oligo", header)
        self.assertNotIn("Reverse Oligo", header)
        self.assertIn("amiRNA*", header)
        # Still one row per result, and no column slipped out of line.
        self.assertEqual(len(rows), 1)
        self.assertEqual(len(rows[0].split("\t")), len(header))

    def test_the_tsv_keeps_the_columns_when_a_vector_is_named(self):
        with tempfile.TemporaryDirectory() as tmp:
            tsv = Path(tmp) / "out.tsv"
            data = _one_result_json(tmp, "out.json",
                                    vector="pMDC32B-AtMIR390a-B/c")
            write_amirna_tsv(data, tsv)
            lines = tsv.read_text().splitlines()
        header = lines[0].split("\t")
        self.assertIn("Forward Oligo", header)
        self.assertEqual(lines[1].split("\t")[header.index("Forward Oligo")],
                         EUDICOT_FWD)

    def test_applying_a_vector_to_a_cached_run_adds_them_back(self):
        """clone_vector.py's path: the cache has no oligos, the result does."""
        with tempfile.TemporaryDirectory() as tmp:
            data = _one_result_json(tmp, "cache.json")
            self.assertNotIn("Forward Oligo", data["results"]["amiRNA Result 1"])
            apply_vector_to_amirna_output(data, "pMDC32B-BS-AtMIR390a-A18G-B/c")
        entry = data["results"]["amiRNA Result 1"]
        self.assertEqual(
            entry["Forward Oligo"],
            make_amirna_oligos(GUIDE, EUDICOT_STAR,
                               "pMDC32B-BS-AtMIR390a-A18G-B/c")[0])
        # The 58-nt pair Alberto reported, rebuilt from a cache holding none.
        self.assertEqual(len(entry["Forward Oligo"]), 58)


class VectorMenuFollowsTheFoldbackTest(unittest.TestCase):
    """
    The menu used to offer all four amiRNA vectors whatever -t said, so
    option 2 (monocot) under the default foldback was a choice that could
    only end in check_vector_foldback's error.
    """

    def _menu(self, foldback, choice="1"):
        printed = []
        with unittest.mock.patch("builtins.input", return_value=choice), \
             unittest.mock.patch("builtins.print", side_effect=printed.append):
            picked = select_vector("amiRNA", foldback)
        return picked, "\n".join(str(p) for p in printed)

    def test_monocot_is_offered_only_its_own_vector(self):
        picked, menu = self._menu("monocot")
        self.assertEqual(picked, "pMDC32B-OsMIR390-B/c")
        for eudicot_vector in (v for v, spec in AMIRNA_VECTORS.items()
                               if spec[3] == "eudicot"):
            self.assertNotIn(eudicot_vector, menu)

    def test_eudicot_is_not_offered_the_monocot_vector(self):
        _, menu = self._menu("eudicot")
        self.assertNotIn("pMDC32B-OsMIR390-B/c", menu)
        self.assertIn("pMDC32B-BS-AtMIR390a-A18G-B/c", menu)

    def test_every_offered_vector_passes_the_check(self):
        """The menu and the gate agree: nothing offered can be refused."""
        for foldback in ("eudicot", "monocot"):
            for choice in range(1, 1 + sum(1 for spec in AMIRNA_VECTORS.values()
                                           if spec[3] == foldback)):
                picked, _ = self._menu(foldback, str(choice))
                check_vector_foldback(picked, foldback)   # must not sys.exit

    def test_no_foldback_still_offers_everything(self):
        """clone_vector.py and older callers pass none."""
        _, menu = self._menu(None)
        for vector in AMIRNA_VECTORS:
            self.assertIn(vector, menu)


class StaleVectorFilesAreFlaggedTest(unittest.TestCase):
    """
    Each vector writes its own file and nothing rewrites the others, so a
    rerun that finds more sites leaves the earlier vector's file holding
    fewer. They are never deleted or rewritten — they are a real result for
    their own vector — but the run says they are there.
    """

    def _note(self, names, vector_suffix="pMDC32B-OsMIR390-Bc"):
        printed = []
        with tempfile.TemporaryDirectory() as tmp:
            for name in names:
                (Path(tmp) / name).write_text("{}")
            with unittest.mock.patch("builtins.print", side_effect=printed.append):
                _warn_stale_vector_files(tmp, "job_amiRNA", vector_suffix)
        return "\n".join(str(p) for p in printed)

    def test_another_vectors_file_is_named(self):
        note = self._note(["job_amiRNA_pMDC32B-AtMIR390a-Bc_psams.json",
                           "job_amiRNA_pMDC32B-OsMIR390-Bc_psams.json"])
        self.assertIn("job_amiRNA_pMDC32B-AtMIR390a-Bc_psams.json", note)
        self.assertNotIn("job_amiRNA_pMDC32B-OsMIR390-Bc_psams.json", note)

    def test_nothing_is_said_when_this_run_is_the_only_vector(self):
        self.assertEqual(
            self._note(["job_amiRNA_pMDC32B-OsMIR390-Bc_psams.json"]), "")

    def test_the_vector_less_copy_is_not_mistaken_for_a_vector(self):
        """{run_key}_psams.json is refreshed in place, not flagged."""
        self.assertEqual(self._note(["job_amiRNA_psams.json"]), "")

    def test_another_runs_files_are_left_alone(self):
        self.assertEqual(
            self._note(["other_amiRNA_pMDC32B-AtMIR390a-Bc_psams.json"]), "")


if __name__ == "__main__":
    unittest.main(verbosity=2)
