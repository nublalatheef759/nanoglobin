from pathlib import Path
import json
import tempfile
import unittest

from nanoglobin.assay import (
    AssayProfileError,
    assay_profile_from_mapping,
    compile_catalogue,
    coverage_rows,
    indistinguishability_classes,
    load_assay_profile,
    read_fasta,
    reverse_complement,
)


FIXTURES = Path(__file__).parent / "fixtures"


class AssayCompilerTests(unittest.TestCase):
    def test_reverse_complement_iupac(self):
        self.assertEqual(reverse_complement("ACGTRYN"), "NRYACGT")

    def test_profile_rejects_unknown_primer(self):
        with self.assertRaises(AssayProfileError):
            assay_profile_from_mapping(
                {
                    "assay_id": "bad",
                    "version": "1",
                    "platform": "ont",
                    "primers": {"F": {"sequence": "ACGT"}},
                    "products": [
                        {
                            "id": "p",
                            "forward_primer": "F",
                            "reverse_primer": "missing",
                            "min_length": 1,
                            "max_length": 100,
                        }
                    ],
                }
            )

    def test_compile_coverage_and_dropout(self):
        profile = load_assay_profile(FIXTURES / "assay.yml")
        haplotypes = read_fasta(FIXTURES / "haplotypes.fa")
        compiled = compile_catalogue(profile, haplotypes)
        by_haplotype = {}
        for product in compiled:
            by_haplotype.setdefault(product.haplotype_id, []).append(product)

        self.assertEqual(len(by_haplotype["alpha_a"]), 1)
        self.assertEqual(len(by_haplotype["alpha_b"]), 1)
        self.assertEqual(len(by_haplotype["alpha_short"]), 1)
        self.assertNotIn("alpha_dropout", by_haplotype)
        self.assertLess(
            by_haplotype["alpha_short"][0].length_bp,
            by_haplotype["alpha_a"][0].length_bp,
        )

        rows = list(coverage_rows(profile, haplotypes, compiled))
        dropout = next(row for row in rows if row["haplotype_id"] == "alpha_dropout")
        self.assertFalse(dropout["observable"])
        self.assertTrue(dropout["required"])

    def test_observationally_identical_haplotypes_are_grouped(self):
        profile = load_assay_profile(FIXTURES / "assay.yml")
        haplotypes = read_fasta(FIXTURES / "haplotypes.fa")
        compiled = compile_catalogue(profile, haplotypes)
        classes = indistinguishability_classes(haplotypes, compiled)
        pair_sets = [set(item.genotype_pairs) for item in classes]
        expected = {
            ("alpha_a", "alpha_a"),
            ("alpha_a", "alpha_b"),
            ("alpha_b", "alpha_b"),
        }
        self.assertTrue(any(expected <= pairs for pairs in pair_sets))

    def test_cli_writes_reproducible_outputs(self):
        from scripts.compile_assay import main

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exit_code = main(
                [
                    "--profile",
                    str(FIXTURES / "assay.yml"),
                    "--haplotypes",
                    str(FIXTURES / "haplotypes.fa"),
                    "--compiled-json",
                    str(root / "compiled.json"),
                    "--coverage-tsv",
                    str(root / "coverage.tsv"),
                    "--indistinguishability-tsv",
                    str(root / "ambiguity.tsv"),
                    "--summary-json",
                    str(root / "summary.json"),
                ]
            )
            self.assertEqual(exit_code, 0)
            payload = json.loads((root / "compiled.json").read_text())
            self.assertEqual(payload["assay"]["assay_id"], "synthetic_globin_panel")
            self.assertTrue((root / "coverage.tsv").is_file())
            self.assertTrue((root / "ambiguity.tsv").is_file())


if __name__ == "__main__":
    unittest.main()
