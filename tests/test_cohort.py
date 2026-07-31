from pathlib import Path
import csv
import tempfile
import unittest

from nanoglobin.cohort import load_cohort_contract, normalise_rows, read_source_csv


FIXTURES = Path(__file__).parent / "fixtures"


class CohortNormalisationTests(unittest.TestCase):
    def test_selective_reporting_is_not_reference(self):
        contract = load_cohort_contract(FIXTURES / "cohort_contract.yml")
        rows = read_source_csv(FIXTURES / "cohort.csv", contract)
        episodes, findings, phenotypes = normalise_rows(
            rows, contract, source_label="fixture"
        )
        self.assertEqual(len(episodes), 2)
        self.assertEqual(len(findings), 4)
        self.assertEqual(len(phenotypes), 3)

        by_id = {finding["finding_id"]: finding for finding in findings}
        self.assertEqual(by_id["E1:HBA"]["reporting_status"], "not_tested")
        self.assertEqual(by_id["E1:HBB"]["reporting_status"], "reported")
        self.assertEqual(by_id["E2:HBA"]["reporting_status"], "not_reported")
        self.assertEqual(by_id["E2:HBB"]["reporting_status"], "not_reported")
        self.assertTrue(
            all(finding["reference_genotype_inferred"] == "false" for finding in findings)
        )
        self.assertEqual(episodes[0]["age_years"], "21")
        self.assertAlmostEqual(float(episodes[1]["age_years"]), 0.5)

    def test_cli_emits_audit_tables(self):
        from scripts.normalise_cohort_reports import main

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exit_code = main(
                [
                    "--input",
                    str(FIXTURES / "cohort.csv"),
                    "--contract",
                    str(FIXTURES / "cohort_contract.yml"),
                    "--episodes",
                    str(root / "episodes.csv"),
                    "--findings",
                    str(root / "findings.csv"),
                    "--phenotypes",
                    str(root / "phenotypes.csv"),
                    "--summary-json",
                    str(root / "summary.json"),
                ]
            )
            self.assertEqual(exit_code, 0)
            with (root / "findings.csv").open(encoding="utf-8-sig") as handle:
                findings = list(csv.DictReader(handle))
            self.assertEqual(len(findings), 4)
            self.assertTrue((root / "summary.json").is_file())


if __name__ == "__main__":
    unittest.main()
