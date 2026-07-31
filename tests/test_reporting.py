from pathlib import Path
import csv
import tempfile
import unittest

from scripts.comprehensive_report import build_report
from scripts.detect_cnv import write_cnv_calls
from scripts.sample_report import main as sample_report_main


class ReportingTests(unittest.TestCase):
    def _write_sample_inputs(self, root: Path, sample: str) -> None:
        results = root / "results"
        variants = root / "variants" / sample
        results.mkdir(parents=True, exist_ok=True)
        variants.mkdir(parents=True, exist_ok=True)
        with (results / f"{sample}.annotated.csv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [
                    "Sample",
                    "Chromosome",
                    "Position",
                    "Ref",
                    "Alt",
                    "HGVS",
                    "Gene",
                    "Consequence",
                    "Genotype",
                    "Quality",
                ]
            )
            writer.writerow(
                [sample, "chr11", 5227002, "A", "T", "HBB:c.20A>T", "HBB", "missense_variant", "0/1", 30]
            )
        vcf_header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        for tool in ("sniffles", "cutesv"):
            (variants / f"{sample}.{tool}.vcf").write_text(vcf_header, encoding="utf-8")
        (variants / f"{sample}.coverage.tsv").write_text(
            "HBA\tmedian_ratio=1.0\n", encoding="utf-8"
        )

    def test_report_contains_only_declared_samples(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write_sample_inputs(root, "REAL")
            self._write_sample_inputs(root, "SIMULATED")
            output = root / "report.csv"
            build_report(
                output,
                ["REAL"],
                results_dir=root / "results",
                variants_dir=root / "variants",
            )
            text = output.read_text(encoding="utf-8-sig")
            self.assertIn("REAL", text)
            self.assertNotIn("SIMULATED", text)

    def test_report_rejects_sample_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write_sample_inputs(root, "REAL")
            annotated = root / "results" / "REAL.annotated.csv"
            text = annotated.read_text()
            annotated.write_text(text.replace("REAL,chr11", "WRONG,chr11"))
            with self.assertRaises(ValueError):
                build_report(
                    root / "report.csv",
                    ["REAL"],
                    results_dir=root / "results",
                    variants_dir=root / "variants",
                )

    def test_coverage_gain_is_candidate_not_pathogenic(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bins = root / "S1.coverage_bins.tsv"
            with bins.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["region", "bin_start", "bin_end", "ratio_to_HBB"])
                for start, ratio in [
                    (170000, 1.0),
                    (170200, 1.7),
                    (170400, 1.8),
                    (170600, 1.7),
                    (170800, 1.0),
                ]:
                    writer.writerow(["HBA", start, start + 199, ratio])
            cnv = root / "cnv.csv"
            write_cnv_calls(cnv, [("S1", bins)])
            with cnv.open(encoding="utf-8-sig") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(rows[0]["Type"], "gain")
            self.assertIn(rows[0]["Naming_status"], {"uncatalogued_gain_candidate", "catalogue_candidate_only"})

            report = root / "report.csv"
            with report.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    [
                        "Sample",
                        "Tool",
                        "Chromosome",
                        "Position",
                        "Ref",
                        "Alt",
                        "HGVS",
                        "Gene",
                        "Consequence",
                        "Genotype",
                        "Quality",
                    ]
                )
                writer.writerow(["S1", "Coverage", "", "", "", "", "", "HBA", "normal", "", ""])
            variants = root / "variants.csv"
            variants.write_text(
                "HVGS,ClinVar classification,Common name,Functionality\n",
                encoding="utf-8",
            )
            sample_report = root / "sample_report.csv"
            self.assertEqual(
                sample_report_main([str(report), str(variants), str(sample_report), str(cnv)]),
                0,
            )
            with sample_report.open(encoding="utf-8-sig") as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(row["N_pathogenic"], "0")
            self.assertEqual(row["N_vus"], "1")
            self.assertIn("gain candidate", row["All_variants_ranked"])


if __name__ == "__main__":
    unittest.main()
