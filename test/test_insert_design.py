#!/usr/bin/env python3
"""
Tests for the full-dsDNA-insert vectors (PVX, TRV and the B/c multimodule
Golden Gate entry), ported from Alberto's Oligo Designer Suite HTML.

The reference values come from the HTML itself: its literal flank and
scaffold sequences, its stated module and insert lengths, and its own
worked example amiR-NbSu.

Run from the repository root with:  python3 test/test_insert_design.py
"""

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.insert_design import (  # noqa: E402
    AMI_COMP,
    AT_MIR173A,
    INSERT_VECTORS,
    ami_star,
    amirna_module,
    build_csv,
    build_fasta,
    build_insert,
    build_txt,
    check_foldback,
    clean_sequence,
    design_basename,
    mir173a_module,
    parse_modules,
    syntasirna_module,
)
from src.utils import oligo_designer  # noqa: E402

# The HTML's worked example, and a second guide from the same preset.
AMI_NBSU = "TGTATGACTCCCGGAATTCCA"
NB_DXS1 = "TAAACCGCGGGTTCCTAACAG"
NB_MIR482A_TS = "GTGGTATGGGGGGAGTCGGGAA"


def star_of(amirna: str) -> str:
    """The amiRNA* exactly as the rest of P-SAMS computes it."""
    star, _fwd, _rev = oligo_designer(amirna, "eudicot")
    return star


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class FlankTest(unittest.TestCase):
    """The flanks must be the HTML's, to the base."""

    def test_flanks_are_the_html_sequences(self):
        self.assertEqual(
            INSERT_VECTORS["PVX"],
            ("AGAGGTCAGCACCAGCTAGC", "AGGGTTTGTTAAGTTTCCCT"),
        )
        self.assertEqual(
            INSERT_VECTORS["TRV"],
            ("CACTTACCCGAGTTAACGCC", "ATGTCCCGAAGACATTAAAC"),
        )
        self.assertEqual(
            INSERT_VECTORS["pMDC32B-B/c-multimodule"],
            ("GCATCGGGTCTCATGTA", "CATTAGAGACCCGATGC"),
        )

    def test_multimodule_flanks_carry_one_bsai_site_per_strand(self):
        design = build_insert(
            "pMDC32B-B/c-multimodule",
            [amirna_module("m", AMI_NBSU, star_of(AMI_NBSU))],
            "bsai",
        )
        # GGTCTC is BsaI; GAGACC is the same site read on the other strand.
        self.assertEqual(design["insert"].count("GGTCTC"), 1)
        self.assertEqual(design["insert"].count("GAGACC"), 1)

    def test_viral_vectors_carry_no_bsai_site(self):
        for vector in ("PVX", "TRV"):
            design = build_insert(
                vector, [amirna_module("m", AMI_NBSU, star_of(AMI_NBSU))], "viral"
            )
            self.assertNotIn("GGTCTC", design["insert"], vector)
            self.assertNotIn("GAGACC", design["insert"], vector)


class ModuleTest(unittest.TestCase):
    """Module geometry, against the lengths the HTML states."""

    def test_amirna_module_is_89_nt_and_pvx_insert_is_129_bp(self):
        module = amirna_module("amiR-NbSu", AMI_NBSU, star_of(AMI_NBSU))
        self.assertEqual(len(module["sequence"]), 89)
        design = build_insert("PVX", [module], "one-amirna")
        self.assertEqual(design["insert_length"], 129)

    def test_amirna_module_matches_the_trv_frames_own_convention(self):
        """
        The TRV app in the HTML splits the duplex differently: it folds two
        scaffold nucleotides into the amiRNA* it reports and uses a 17-nt 3'
        basal stem instead of a 15-nt one. The DNA that comes out is the
        same, and this pins that down, since P-SAMS reports the other
        convention and must not silently produce a different construct.
        """
        special = {"A": "C", "G": "T", "C": "A", "T": "G"}
        rev = revcomp(AMI_NBSU)
        trv_star = rev[:10] + special[rev[10]] + rev[11:21]
        trv_module = (
            "AGTAGAGAAGAATCTGTA" + AMI_NBSU + "CGAAATCAAACT"
            + trv_star + "CATTGGCTCTTCTTACT"
        )
        ours = amirna_module("m", AMI_NBSU, star_of(AMI_NBSU))["sequence"]
        self.assertEqual(ours, trv_module)

    def test_syntasirna_module_is_target_site_spacer_then_guides(self):
        module = syntasirna_module(
            "NbmiR482a", NB_MIR482A_TS, [("a", AMI_NBSU), ("b", NB_DXS1)]
        )
        self.assertEqual(
            module["sequence"], NB_MIR482A_TS + "TAGACCATTTA" + AMI_NBSU + NB_DXS1
        )
        self.assertEqual(len(module["sequence"]), 22 + 11 + 21 + 21)

    def test_mir173a_module_is_the_fixed_precursor_verbatim(self):
        # Transcribed from the HTML's own AT_MIR173A constant, not retyped.
        self.assertEqual(
            AT_MIR173A,
            "TGGTGATTAAGTACTTTCGCTTGCAGAGAGAAATCACAGTGGTCAAAAAAGTTGTAGTTTT"
            "CTTAAAGTCTCTTTCCTCTGTGATTCTCTGTGTAAGCGAAAGAGCTTGCTCCCTAA",
        )
        self.assertEqual(len(AT_MIR173A), 117)
        self.assertEqual(mir173a_module()["sequence"], AT_MIR173A)

    def test_modules_are_butt_joined_with_no_spacer_between_them(self):
        a = amirna_module("a", AMI_NBSU, star_of(AMI_NBSU))
        b = syntasirna_module("NbmiR482a", NB_MIR482A_TS, [("g", NB_DXS1)])
        design = build_insert("PVX", [a, b], "tandem")
        flank5, flank3 = INSERT_VECTORS["PVX"]
        self.assertEqual(
            design["insert"],
            flank5 + a["sequence"] + b["sequence"] + flank3,
        )

    def test_empty_module_list_is_refused(self):
        with self.assertRaises(SystemExit):
            build_insert("PVX", [], "nothing")


class ModuleSpecTest(unittest.TestCase):
    """Modules are addressed by their cached result, never pasted in."""

    def setUp(self):
        self.amirna_index = {
            1: {"amiRNA": AMI_NBSU, "amiRNA*": star_of(AMI_NBSU)},
            2: {"amiRNA": NB_DXS1, "amiRNA*": star_of(NB_DXS1)},
        }
        self.syn_index = {"1.1": AMI_NBSU, "2.1": NB_DXS1}

    def test_hybrid_spec_resolves_in_the_order_given(self):
        modules = parse_modules(
            "ami:2,syn:NbmiR482a:1.1+2.1", self.amirna_index, self.syn_index
        )
        self.assertEqual([m["type"] for m in modules], ["amiRNA", "syn-tasiRNA"])
        self.assertEqual(modules[0]["amiRNA"], NB_DXS1)
        self.assertEqual(modules[1]["target_site"], NB_MIR482A_TS)
        self.assertEqual(
            [seq for _name, seq in modules[1]["guides"]], [AMI_NBSU, NB_DXS1]
        )

    def test_target_site_may_be_spelled_out_instead_of_named(self):
        modules = parse_modules(
            f"syn:{NB_MIR482A_TS}:1.1", self.amirna_index, self.syn_index
        )
        self.assertEqual(modules[0]["target_site"], NB_MIR482A_TS)

    def test_unknown_references_are_refused(self):
        for spec in ("ami:9", "syn:NbmiR482a:9.9", "syn:nope:1.1", "junk:1"):
            with self.subTest(spec=spec):
                with self.assertRaises(SystemExit):
                    parse_modules(spec, self.amirna_index, self.syn_index)


class ExportTest(unittest.TestCase):
    def setUp(self):
        self.design = build_insert(
            "PVX",
            [
                amirna_module("amiR-NbSu", AMI_NBSU, star_of(AMI_NBSU)),
                syntasirna_module("NbmiR482a", NB_MIR482A_TS, [("g", NB_DXS1)]),
            ],
            "demo construct",
        )

    def test_fasta_is_one_header_and_one_unwrapped_sequence(self):
        lines = build_fasta(self.design).splitlines()
        self.assertEqual(len(lines), 2)
        self.assertEqual(lines[0], ">demo_construct_full_insert")
        self.assertEqual(lines[1], self.design["insert"])

    def test_csv_is_a_field_value_table_carrying_the_insert(self):
        rows = list(csv.reader(build_csv(self.design).splitlines()))
        self.assertEqual(rows[0], ["field", "value"])
        fields = dict(rows[1:])
        self.assertEqual(fields["full_insert_5to3"], self.design["insert"])
        self.assertEqual(fields["module_1_amiRNA"], AMI_NBSU)
        self.assertEqual(fields["module_2_target_seq"], NB_MIR482A_TS)

    def test_txt_reports_the_insert_and_both_flanks(self):
        text = build_txt(self.design)
        self.assertIn(self.design["insert"], text)
        self.assertIn(INSERT_VECTORS["PVX"][0], text)
        self.assertIn(INSERT_VECTORS["PVX"][1], text)
        self.assertIn("PVX hybrid (amiRNA + syn-tasiRNA in tandem)", text)


class CommandLineTest(unittest.TestCase):
    """clone_vector.py builds an insert from a previous run's cache."""

    def _make_run(self, tmp: Path) -> Path:
        cache = tmp / ".cache"
        cache.mkdir(parents=True)
        (cache / "ACC_amiRNA_psams.json").write_text(json.dumps({
            "optimal": {
                "amiRNA Optimal Result 1": {
                    "amiRNA": AMI_NBSU, "amiRNA*": star_of(AMI_NBSU),
                    "Forward Oligo": "", "Reverse Oligo": "", "TargetFinder": [],
                },
            },
            "suboptimal": {},
        }))
        (cache / "ACC_syntasiRNA_psams.json").write_text(json.dumps({
            "blocks": [{
                "name": "Gene set 1",
                "optimal": {"optimal 1.1": {"syn-tasiRNA": NB_DXS1, "TargetFinder": []}},
                "suboptimal": {},
            }],
        }))
        return tmp

    def test_hybrid_insert_from_cached_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self._make_run(Path(tmp))
            result = subprocess.run(
                [sys.executable, "clone_vector.py",
                 "-o", str(run), "-V", "TRV",
                 "-M", "ami:1,syn:NbmiR482a:1.1", "-n", "hybrid"],
                cwd=REPO_ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            base = run / "hybrid_TRV_multimodule"
            for suffix in ("_design.json", "_design.txt",
                           "_design.fasta", "_summary.csv"):
                self.assertTrue(Path(str(base) + suffix).exists(), suffix)

            design = json.loads(Path(str(base) + "_design.json").read_text())
            self.assertEqual(design["vector"], "TRV")
            self.assertEqual([m["type"] for m in design["modules"]],
                             ["amiRNA", "syn-tasiRNA"])
            # 20 + 89 + (22 + 11 + 21) + 20
            self.assertEqual(design["insert_length"], 183)

    def test_modules_and_oligo_pair_vectors_do_not_mix(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self._make_run(Path(tmp))
            result = subprocess.run(
                [sys.executable, "clone_vector.py",
                 "-o", str(run), "-V", "pMDC32B-AtMIR390a-B/c", "-M", "ami:1"],
                cwd=REPO_ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("--modules only applies", result.stderr)


class StandaloneSequenceTest(unittest.TestCase):
    """Sequences typed straight into the command line, with no run behind them."""

    def test_derived_star_matches_the_html_worked_example(self):
        self.assertEqual(ami_star(AMI_NBSU), "GAATTCCGTGAGTCATACACA")

    def test_derived_star_matches_the_trv_frames_own_convention(self):
        """
        The TRV app derives the star as the plain reverse complement with
        position 10 swapped through a different table. It must never disagree
        with the convention used here.
        """
        special = {"A": "C", "G": "T", "C": "A", "T": "G"}
        for amirna in (AMI_NBSU, NB_DXS1, "ACGTACGTACGTACGTACGTA"):
            rc = "".join(AMI_COMP[b] for b in reversed(amirna))
            trv = rc[:10] + special[rc[10]] + rc[11:]
            # The TRV frame puts the 2-nt scaffold dinucleotide at the front
            # of its star and the trailing "CA" in the basal stem, so the 19 nt
            # the two conventions share start at index 2 there.
            self.assertEqual(ami_star(amirna)[:19], trv[2:21], amirna)

    def test_derived_star_is_reused_when_one_is_not_given(self):
        typed = parse_modules(f"ami:x:{AMI_NBSU}")[0]
        spelled = parse_modules(f"ami:x:{AMI_NBSU}:{ami_star(AMI_NBSU)}")[0]
        self.assertEqual(typed["sequence"], spelled["sequence"])

    def test_sequences_are_cleaned_the_way_the_html_cleans_them(self):
        self.assertEqual(
            clean_sequence(" ugua-tgactcccggaattcca ", "amiRNA"), AMI_NBSU)

    def test_non_acgt_input_is_refused(self):
        with self.assertRaises(SystemExit):
            clean_sequence("ACGTX", "amiRNA")

    def test_wrong_lengths_are_refused(self):
        with self.assertRaises(SystemExit):
            amirna_module("x", "ACGT")
        with self.assertRaises(SystemExit):
            syntasirna_module("ts", "ACGT", [("g", NB_DXS1)])
        with self.assertRaises(SystemExit):
            syntasirna_module("ts", NB_MIR482A_TS, [("g", "ACGT")])

    def test_a_module_needs_at_least_one_syn_tasirna(self):
        with self.assertRaises(SystemExit):
            syntasirna_module("ts", NB_MIR482A_TS, [])

    def test_cached_and_typed_out_guides_mix_in_one_module(self):
        module = parse_modules(
            f"syn:NbmiR482a:1.1+NbDXS1={NB_DXS1}",
            {}, {"1.1": AMI_NBSU},
        )[0]
        self.assertEqual([name for name, _ in module["guides"]],
                         ["syn-tasiRNA 1.1", "NbDXS1"])
        self.assertEqual([seq for _, seq in module["guides"]],
                         [AMI_NBSU, NB_DXS1])

    def test_target_site_may_carry_its_own_name(self):
        module = parse_modules(f"syn:MyTS={NB_MIR482A_TS}:g={NB_DXS1}")[0]
        self.assertEqual(module["target_name"], "MyTS")
        self.assertEqual(module["target_site"], NB_MIR482A_TS)

    def test_mir173a_module_may_be_renamed(self):
        self.assertEqual(parse_modules("mir173a:my173")[0]["name"], "my173")


class HtmlParityTest(unittest.TestCase):
    """The labels, wording and filenames the HTML produces, reproduced here."""

    def _design(self, vector, spec, name):
        return build_insert(vector, parse_modules(spec), name)

    def test_type_labels_match_the_html(self):
        ami = f"ami:a:{AMI_NBSU}"
        syn = f"syn:NbmiR482a:g={NB_DXS1}"
        cases = [
            ("PVX", ami, "PVX amiRNA (shc / AtMIR390a-based)"),
            ("PVX", syn, "PVX syn-tasiRNA (AtTAS1c-based)"),
            ("PVX", f"{ami},{syn}", "PVX hybrid (amiRNA + syn-tasiRNA in tandem)"),
            ("TRV", ami, "TRV amiRNA (shc / MIR390-based)"),
            ("TRV", syn, "TRV syn-tasiRNA (AtTAS1c-based)"),
            ("TRV", f"{ami},{syn}", "TRV Multimodule (amiRNA + syn-tasiRNA in tandem)"),
            ("pMDC32B-B/c-multimodule", f"{ami},{syn}",
             "B/c multimodule (amiRNA + syn-tasiRNA in tandem)"),
        ]
        for vector, spec, expected in cases:
            self.assertEqual(self._design(vector, spec, "x")["type"], expected, spec)

    def test_an_atmir173a_only_insert_is_not_called_a_mir390_amirna(self):
        design = self._design("PVX", "mir173a", "x")
        self.assertNotIn("MIR390", design["type"])

    def test_the_bc_insert_is_reported_as_single_stranded(self):
        ami = f"ami:a:{AMI_NBSU}"
        self.assertIn("Full ssDNA insert",
                      build_txt(self._design("pMDC32B-B/c-multimodule", ami, "x")))
        self.assertIn("Full dsDNA insert",
                      build_txt(self._design("PVX", ami, "x")))

    def test_built_in_target_sites_carry_the_html_label(self):
        module = parse_modules(f"syn:NbmiR482a:g={NB_DXS1}")[0]
        self.assertEqual(module["target_name"], "NbmiR482a TS")

    def test_filenames_follow_the_html_downloads(self):
        design = self._design("TRV", f"ami:a:{AMI_NBSU}", "TRV_amiRNA_NbSu")
        self.assertEqual(design_basename(design), "TRV_amiRNA_NbSu_TRV_amiRNA")

    def test_an_unnamed_construct_falls_back_the_way_the_html_does(self):
        design = self._design("PVX", f"syn:NbmiR482a:g={NB_DXS1}", "")
        self.assertEqual(design_basename(design),
                         "pvx_syn_tasiRNA_PVX_syn_tasiRNA")

    def test_worked_example_lengths_match_the_html(self):
        ami = f"ami:amiR-NbSu:{AMI_NBSU}"
        syn = f"syn:NbmiR482a:amiR-NbSu={AMI_NBSU}+NbDXS1={NB_DXS1}"
        self.assertEqual(self._design("PVX", ami, "x")["insert_length"], 129)
        self.assertEqual(self._design("TRV", ami, "x")["insert_length"], 129)
        self.assertEqual(self._design("PVX", syn, "x")["insert_length"], 115)
        self.assertEqual(self._design("PVX", f"{ami},{syn}", "x")["insert_length"], 204)
        self.assertEqual(
            self._design("pMDC32B-B/c-multimodule", f"{ami},{syn}", "x")["insert_length"],
            198)

    def test_bc_worked_example_is_byte_identical_to_the_html(self):
        design = self._design(
            "pMDC32B-B/c-multimodule",
            f"ami:amiR-NbSu:{AMI_NBSU},"
            f"syn:NbmiR482a:amiR-NbSu={AMI_NBSU}+NbDXS1={NB_DXS1}",
            "Bc_multimodule_example")
        self.assertEqual(design["insert"], (
            "GCATCGGGTCTCATGTAAGTAGAGAAGAATCTGTATGTATGACTCCCGGAATTCCACGAAAT"
            "CAAACTTGGAATTCCGTGAGTCATACACATTGGCTCTTCTTACTGTGGTATGGGGGGAGTCG"
            "GGAATAGACCATTTATGTATGACTCCCGGAATTCCATAAACCGCGGGTTCCTAACAGCATTA"
            "GAGACCCGATGC"))


class FoldbackGuardTest(unittest.TestCase):
    """Only the eudicot AtMIR390a foldback has a defined insert design."""

    def test_a_monocot_run_is_refused(self):
        with self.assertRaises(SystemExit):
            check_foldback("monocot")

    def test_eudicot_and_pre_existing_caches_are_accepted(self):
        check_foldback("eudicot")
        check_foldback(None)


class StandaloneCommandLineTest(unittest.TestCase):
    """clone_vector.py with no cached run anywhere in sight."""

    def test_insert_is_built_without_any_cached_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = subprocess.run(
                [sys.executable, "clone_vector.py",
                 "-o", tmp, "-V", "PVX",
                 "-M", f"ami:amiR-NbSu:{AMI_NBSU}", "-n", "PVX_amiRNA_NbSu"],
                cwd=REPO_ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("insert 129 bp", result.stdout)

            base = Path(tmp) / "PVX_amiRNA_NbSu_PVX_amiRNA"
            design = json.loads(Path(str(base) + "_design.json").read_text())
            self.assertEqual(design["modules"][0]["amiRNA*"], ami_star(AMI_NBSU))
            self.assertIn("Full dsDNA insert",
                          Path(str(base) + "_design.txt").read_text())

    def test_a_monocot_cache_is_refused_on_the_command_line(self):
        with tempfile.TemporaryDirectory() as tmp:
            cache = Path(tmp) / ".cache"
            cache.mkdir()
            (cache / "ACC_amiRNA_psams.json").write_text(json.dumps({
                "foldback": "monocot",
                "optimal": {"amiRNA Optimal Result 1": {
                    "amiRNA": AMI_NBSU, "amiRNA*": star_of(AMI_NBSU),
                    "Forward Oligo": "", "Reverse Oligo": "", "TargetFinder": [],
                }},
                "suboptimal": {},
            }))
            result = subprocess.run(
                [sys.executable, "clone_vector.py",
                 "-o", tmp, "-V", "PVX", "-M", "ami:1"],
                cwd=REPO_ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("monocot", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
