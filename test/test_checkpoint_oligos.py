#!/usr/bin/env python3
"""
Regression test for the two disagreeing copies of the amiRNA cloning oligos.

A finished amiRNA job left two oligo tables in the same output folder:

  .cache/<run_key>_optimal_results.tsv   Oligo1/Oligo2        75 nt
  <run_key>_<vector>_psams.tsv           Forward/Reverse      58 nt

The second one is the one to order. The first held the bare foldback
hairpin from oligo_designer(), written before any vector is applied, and
nothing in its name said so — open the checkpoint and you order the wrong
construct. The checkpoint now carries no oligos at all, and a resumed run
recomputes them from the guide instead of reading them back.

Run from the repository root with:  python3 test/test_checkpoint_oligos.py
"""

import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.oligo_design import AMIRNA_VECTORS, make_amirna_oligos  # noqa: E402
from src.pipeline import serial_jobs  # noqa: E402
from src.utils import Background, create_outputs, load_resume_state, oligo_designer  # noqa: E402
# The same throwaway transcript/guide fixtures the pepper tests build, so a
# checkpoint row can be produced without real sequence data.
from test_pepper_bugs import (  # noqa: E402
    GUIDE, TARGETFINDER, TARGET_SITE, TRANSCRIPT, make_annotation_db, write_db)

FOLDBACK = "eudicot"


class CheckpointHasNoOligoColumnsTest(unittest.TestCase):

    def _headers(self, syntasirna):
        with tempfile.TemporaryDirectory() as tmp:
            subopt_name, opt_name = create_outputs("job_amiRNA", tmp,
                                                   syntasirna=syntasirna)
            return (Path(opt_name).read_text().rstrip("\n").split("\t"),
                    Path(subopt_name).read_text().rstrip("\n").split("\t"))

    def test_amirna_checkpoint_has_no_oligo_columns(self):
        opt, subopt = self._headers(syntasirna=False)
        self.assertEqual(
            opt, ["Site_index", "Guide", "Star", "Isoforms", "Isoform_seqs"])
        self.assertEqual(
            subopt, ["Site_index", "Offtarget_N", "Offtarget_list", "Guide",
                     "Star", "Isoforms", "Isoform_seqs"])

    def test_syntasirna_checkpoint_keeps_its_gene_set_column(self):
        opt, subopt = self._headers(syntasirna=True)
        self.assertEqual(opt[0], "Gene_set")
        self.assertEqual(subopt[0], "Gene_set")
        for header in (opt, subopt):
            self.assertNotIn("Oligo1", header)
            self.assertNotIn("Oligo2", header)

    def test_a_row_the_pipeline_writes_matches_its_header(self):
        """
        The header and the row builders live in different modules, so a
        column dropped from one and not the other would go unnoticed until
        a resumed run read a guide out of the Star column. Run the real
        pipeline on a throwaway transcript and compare the row it appends
        against the header create_outputs wrote.
        """
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            # Two identical transcripts: the second is an unavoidable
            # off-target, so the candidate lands in the suboptimal file.
            db = write_db(tmp / "db.fa", ["TX00001", "TX00002"])
            annot = make_annotation_db(
                tmp / "annot.db",
                [("TX00001", "intended target"), ("TX00002", "off-target")])
            conn = sqlite3.connect(annot)
            conn.row_factory = sqlite3.Row
            tf_dir = tmp / "tf_results"
            tf_dir.mkdir()

            subopt_name, opt_name = create_outputs("job_amiRNA", tmp)
            site = {
                "guide": GUIDE, "star": GUIDE,
                "oligo1": "SHOULD-NOT-BE-WRITTEN",
                "oligo2": "SHOULD-NOT-BE-WRITTEN",
                "names": "TX00001", "seqs": TARGET_SITE,
            }
            serial_jobs(
                target_count=1, construct="amiRNA",
                ids={"TX00001": TRANSCRIPT}, site_scores=[site],
                targetfinder=str(TARGETFINDER),
                backgrounds=[Background("test_species", conn, str(db))],
                bg=True, subopt_name=subopt_name, opt_name=opt_name,
                output_folder=str(tmp), accession_list=["TX00001"],
                tf_dir=str(tf_dir), potential_target_n=1, fasta=True)
            conn.close()

            lines = Path(subopt_name).read_text().splitlines()

        self.assertEqual(len(lines), 2, "expected one suboptimal row")
        header, row = (line.split("\t") for line in lines)
        self.assertEqual(len(row), len(header))
        self.assertNotIn("SHOULD-NOT-BE-WRITTEN", "\t".join(row))
        self.assertEqual(row[header.index("Guide")], GUIDE)


class ResumeRecomputesOligosTest(unittest.TestCase):

    OPT_HEADER = "Site_index\tGuide\tStar\tIsoforms\tIsoform_seqs\n"
    SUB_HEADER = ("Site_index\tOfftarget_N\tOfftarget_list\tGuide\tStar\t"
                  "Isoforms\tIsoform_seqs\n")

    def _resume(self, tmp, **kwargs):
        tf_dir = Path(tmp) / "tf_results"
        tf_dir.mkdir()
        (tf_dir / "site_0001_TargetFinder_result.json").write_text('{"ok": 1}')
        (tf_dir / "site_0002_TargetFinder_result.json").write_text('{"ok": 2}')

        star, _, _ = oligo_designer(GUIDE, FOLDBACK)
        opt = Path(tmp) / "opt.tsv"
        sub = Path(tmp) / "sub.tsv"
        opt.write_text(self.OPT_HEADER + f"1\t{GUIDE}\t{star}\tn\ts\n")
        sub.write_text(self.SUB_HEADER + f"2\t3\tx,y\t{GUIDE}\t{star}\tn\ts\n")

        return load_resume_state(str(tf_dir), str(opt), str(sub), **kwargs)

    def test_amirna_resume_rebuilds_the_foldback_pair_from_the_guide(self):
        expected_star, expected1, expected2 = oligo_designer(GUIDE, FOLDBACK)
        with tempfile.TemporaryDirectory() as tmp:
            opt, subopt, seen, start = self._resume(tmp, foldback=FOLDBACK)

        self.assertEqual(opt[0]["oligo1"], expected1)
        self.assertEqual(opt[0]["oligo2"], expected2)
        self.assertEqual(subopt[0]["site"]["oligo1"], expected1)
        self.assertEqual(subopt[0]["site"]["oligo2"], expected2)
        # Everything the resume already did keeps working.
        self.assertEqual(opt[0]["star"], expected_star)
        self.assertEqual(seen, {GUIDE})
        self.assertEqual(start, 2)

    def test_syntasirna_resume_leaves_the_oligos_empty(self):
        """No foldback is passed for syntasiRNA: it has no per-site pair."""
        with tempfile.TemporaryDirectory() as tmp:
            opt, subopt, _, _ = self._resume(tmp)

        self.assertEqual(opt[0]["oligo1"], "")
        self.assertEqual(subopt[0]["site"]["oligo2"], "")

    def test_an_old_checkpoint_with_oligo_columns_still_resumes(self):
        """
        Folders written by an earlier version still carry Oligo1/Oligo2.
        Rows are read by header name, so the extra columns are ignored and
        the oligos are recomputed — the stale 75-nt pair is never reused.
        """
        star, expected1, _ = oligo_designer(GUIDE, FOLDBACK)
        with tempfile.TemporaryDirectory() as tmp:
            tf_dir = Path(tmp) / "tf_results"
            tf_dir.mkdir()
            (tf_dir / "site_0001_TargetFinder_result.json").write_text('{"ok": 1}')
            opt = Path(tmp) / "opt.tsv"
            sub = Path(tmp) / "sub.tsv"
            opt.write_text(
                "Site_index\tGuide\tStar\tOligo1\tOligo2\tIsoforms\tIsoform_seqs\n"
                f"1\t{GUIDE}\t{star}\tSTALE1\tSTALE2\tn\ts\n")
            sub.write_text(self.SUB_HEADER)

            opt_sites, _, seen, _ = load_resume_state(
                str(tf_dir), str(opt), str(sub), foldback=FOLDBACK)

        self.assertEqual(opt_sites[0]["guide"], GUIDE)
        self.assertEqual(opt_sites[0]["oligo1"], expected1)
        self.assertEqual(seen, {GUIDE})


class VectorOligosAreTheOnlyOnesTest(unittest.TestCase):
    """
    Why no pair written before the vector is known can be right: picking a
    vector changes the BsaI prefix and the spacer, so it changes both the
    sequence and the length of the oligos.
    """

    def test_only_the_plain_atmir390a_vector_matches_the_bare_foldback_pair(self):
        """
        pMDC32B-AtMIR390a-B/c reuses the eudicot foldback verbatim (same
        TGTA prefix, same 29-nt loop), so for that one vector the two
        tables happened to agree — which is how the mismatch stayed
        invisible until a BS-AtMIR390a job was run. Every other vector
        disagrees.
        """
        star, foldback1, _ = oligo_designer(GUIDE, FOLDBACK)
        agreeing = [v for v in AMIRNA_VECTORS
                    if make_amirna_oligos(GUIDE, star, v)[0] == foldback1]
        self.assertEqual(agreeing, ["pMDC32B-AtMIR390a-B/c"])

    def test_the_a18g_vector_is_the_58_nt_pair_alberto_reported(self):
        star, foldback1, _ = oligo_designer(GUIDE, FOLDBACK)
        fwd, rev = make_amirna_oligos(
            GUIDE, star, "pMDC32B-BS-AtMIR390a-A18G-B/c")
        self.assertEqual(len(foldback1), 75)
        self.assertEqual(len(fwd), 58)
        self.assertEqual(len(rev), 58)


if __name__ == "__main__":
    unittest.main(verbosity=2)
