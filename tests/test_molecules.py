from pathlib import Path
import csv
import gzip
import json
import tempfile
import unittest

from nanoglobin.assay import (
    Haplotype,
    assay_profile_from_mapping,
    compile_catalogue,
    load_assay_profile,
    read_fasta,
    reverse_complement,
    write_compiled_json,
)
from nanoglobin.molecules import (
    AdmissionAccumulator,
    AdmissionConfig,
    FastqRecord,
    MoleculeAdmissionError,
    MoleculeObservation,
    admit_molecule,
    load_compiled_assay,
    read_fastq,
)


FIXTURES = Path(__file__).parent / "fixtures"


class MoleculeAdmissionTests(unittest.TestCase):
    def setUp(self):
        self.profile = load_assay_profile(FIXTURES / "assay.yml")
        self.haplotypes = read_fasta(FIXTURES / "haplotypes.fa")
        self.products = compile_catalogue(self.profile, self.haplotypes)
        self.tmp = tempfile.TemporaryDirectory()
        compiled_path = Path(self.tmp.name) / "compiled.json"
        write_compiled_json(
            compiled_path, self.profile, self.haplotypes, self.products
        )
        self.compiled = load_compiled_assay(compiled_path, self.profile)

    def tearDown(self):
        self.tmp.cleanup()

    def test_complete_read_preserves_equivalent_haplotype_ambiguity(self):
        sequence = next(
            product.sequence
            for product in self.products
            if product.haplotype_id == "alpha_a"
        )
        observation = admit_molecule(
            FastqRecord("read-a", sequence, "I" * len(sequence)),
            sample="S1",
            profile=self.profile,
            compiled=self.compiled,
        )
        self.assertEqual(observation.status, "complete")
        self.assertEqual(observation.orientation, "forward")
        self.assertEqual(observation.assigned_product_id, "HBA_product")
        self.assertEqual(observation.assignment_state, "equivalent_haplotypes")
        self.assertEqual(
            observation.candidate_haplotypes, ("alpha_a", "alpha_b")
        )
        self.assertEqual(observation.edit_distance, 0)

    def test_reverse_complement_read_is_oriented(self):
        sequence = next(
            product.sequence
            for product in self.products
            if product.haplotype_id == "alpha_short"
        )
        reverse = reverse_complement(sequence)
        observation = admit_molecule(
            FastqRecord("read-r", reverse, "I" * len(reverse)),
            sample="S1",
            profile=self.profile,
            compiled=self.compiled,
        )
        self.assertEqual(observation.status, "complete")
        self.assertEqual(observation.orientation, "reverse")
        self.assertEqual(observation.assignment_state, "unique_sequence")
        self.assertEqual(observation.candidate_haplotypes, ("alpha_short",))

    def test_one_ended_and_off_target_are_not_reference_calls(self):
        forward = self.profile.primers["HBA_F"].sequence + "AAAAAA"
        one_ended = admit_molecule(
            FastqRecord("one", forward, "I" * len(forward)),
            sample="S1",
            profile=self.profile,
            compiled=self.compiled,
        )
        off_target_sequence = "G" * 24
        off_target = admit_molecule(
            FastqRecord("off", off_target_sequence, "I" * len(off_target_sequence)),
            sample="S1",
            profile=self.profile,
            compiled=self.compiled,
        )
        self.assertEqual(one_ended.status, "one_ended")
        self.assertFalse(one_ended.assigned_product_id)
        self.assertEqual(off_target.status, "off_target")
        self.assertFalse(off_target.assigned_product_id)

    def test_incompatible_terminal_primers_are_chimera_candidate(self):
        profile = assay_profile_from_mapping(
            {
                "assay_id": "two_product_panel",
                "version": "1",
                "platform": "ont",
                "primers": {
                    "AF": {"sequence": "AAAACCCC"},
                    "AR": {"sequence": "GGGGTTTT"},
                    "BF": {"sequence": "ACGTACGT"},
                    "BR": {"sequence": "TGCATGCA"},
                },
                "products": [
                    {
                        "id": "A",
                        "forward_primer": "AF",
                        "reverse_primer": "AR",
                        "min_length": 20,
                        "max_length": 100,
                    },
                    {
                        "id": "B",
                        "forward_primer": "BF",
                        "reverse_primer": "BR",
                        "min_length": 20,
                        "max_length": 100,
                    },
                ],
            }
        )
        haplotypes = [
            Haplotype(
                "ha", "AAAACCCC" + "G" * 20 + reverse_complement("GGGGTTTT")
            ),
            Haplotype(
                "hb", "ACGTACGT" + "C" * 20 + reverse_complement("TGCATGCA")
            ),
        ]
        products = compile_catalogue(profile, haplotypes)
        path = Path(self.tmp.name) / "two.json"
        write_compiled_json(path, profile, haplotypes, products)
        compiled = load_compiled_assay(path, profile)
        chimera = "AAAACCCC" + "T" * 20 + reverse_complement("TGCATGCA")
        observation = admit_molecule(
            FastqRecord("chimera", chimera, "I" * len(chimera)),
            sample="S1",
            profile=profile,
            compiled=compiled,
            config=AdmissionConfig(
                min_primer_overlap=8, max_primer_mismatches=1
            ),
        )
        self.assertEqual(observation.status, "chimera_candidate")
        self.assertEqual(observation.left_primer_id, "AF")
        self.assertEqual(observation.right_primer_id, "BR")

    def test_declared_pair_with_too_short_span_is_primer_dimer_candidate(self):
        forward = self.profile.primers["HBA_F"].sequence
        reverse_binding = reverse_complement(
            self.profile.primers["HBA_R"].sequence
        )
        sequence = forward + reverse_binding
        observation = admit_molecule(
            FastqRecord("dimer", sequence, "I" * len(sequence)),
            sample="S1",
            profile=self.profile,
            compiled=self.compiled,
        )
        self.assertEqual(observation.status, "primer_dimer_candidate")
        self.assertFalse(observation.assigned_product_id)

    def test_primer_seed_index_falls_back_to_exhaustive_matching(self):
        profile = assay_profile_from_mapping(
            {
                "assay_id": "seed_fallback",
                "version": "1",
                "platform": "ont",
                "primers": {
                    "F": {"sequence": "AAAAAAACCCCCCC", "max_mismatches": 0},
                    "R": {"sequence": "TTTTGGGGTTTT", "max_mismatches": 0},
                    "DECOY_F": {"sequence": "GGGGGGGAAAAAAA", "max_mismatches": 0},
                    "DECOY_R": {"sequence": "CCCCCCCTTTTTTT", "max_mismatches": 0},
                },
                "products": [
                    {
                        "id": "true",
                        "forward_primer": "F",
                        "reverse_primer": "R",
                        "min_length": 35,
                        "max_length": 100,
                    },
                    {
                        "id": "decoy",
                        "forward_primer": "DECOY_F",
                        "reverse_primer": "DECOY_R",
                        "min_length": 35,
                        "max_length": 100,
                    },
                ],
            }
        )
        true_template = (
            profile.primers["F"].sequence
            + "GGGGGGGTTTTTTTTTT"
            + reverse_complement(profile.primers["R"].sequence)
        )
        haplotypes = [Haplotype("h", true_template)]
        products = compile_catalogue(profile, haplotypes)
        path = Path(self.tmp.name) / "fallback.json"
        write_compiled_json(path, profile, haplotypes, products)
        compiled = load_compiled_assay(path, profile)
        # Two substitutions disrupt every exact 7-mer in the true forward
        # primer, while an interior decoy seed keeps the shortlist non-empty.
        noisy = list(products[0].sequence)
        noisy[3] = "C"
        noisy[10] = "A"
        noisy = "".join(noisy)
        from nanoglobin.molecules import AdmissionEngine

        engine = AdmissionEngine(
            profile=profile,
            compiled=compiled,
            config=AdmissionConfig(
                primer_seed_length=7,
                min_primer_overlap=10,
                max_primer_mismatches=3,
                max_primer_error_rate=0.30,
            ),
        )
        observation = engine.admit(
            FastqRecord("fallback", noisy, "I" * len(noisy)), sample="S1"
        )
        self.assertEqual(observation.status, "complete")
        self.assertEqual(observation.assigned_product_id, "true")

    def test_accumulator_uses_exact_length_histogram(self):
        accumulator = AdmissionAccumulator()
        for index, length in enumerate((40, 40, 50, 60)):
            accumulator.add(
                MoleculeObservation(
                    sample="S1",
                    read_id=f"r{index}",
                    read_length_bp=length,
                    status="complete",
                    orientation="forward",
                    assigned_product_id="HBA_product",
                    assignment_state="unique_sequence",
                    observed_product_length_bp=length,
                )
            )
        row = next(accumulator.product_rows())
        self.assertEqual(row["median_observed_length_bp"], "45.0")
        self.assertEqual(row["min_observed_length_bp"], 40)
        self.assertEqual(row["max_observed_length_bp"], 60)
        self.assertEqual(accumulator.product_lengths["HBA_product"][40], 2)

    def test_fastq_gzip_and_cli_outputs(self):
        from scripts.admit_amplicon_reads import main

        sequence = next(
            product.sequence
            for product in self.products
            if product.haplotype_id == "alpha_a"
        )
        root = Path(self.tmp.name)
        fastq = root / "reads.fastq.gz"
        with gzip.open(fastq, "wt", encoding="utf-8") as handle:
            handle.write(f"@r1\n{sequence}\n+\n{'I' * len(sequence)}\n")
            handle.write("@off\nGGGGGGGGGGGGGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIIIIIIIIIII\n")
        compiled_path = root / "compiled-cli.json"
        write_compiled_json(
            compiled_path, self.profile, self.haplotypes, self.products
        )
        exit_code = main(
            [
                "--sample",
                "S1",
                "--fastq",
                str(fastq),
                "--profile",
                str(FIXTURES / "assay.yml"),
                "--compiled-json",
                str(compiled_path),
                "--molecules-tsv",
                str(root / "molecules.tsv"),
                "--product-counts-tsv",
                str(root / "products.tsv"),
                "--summary-json",
                str(root / "summary.json"),
            ]
        )
        self.assertEqual(exit_code, 0)
        with (root / "molecules.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["status"] for row in rows], ["complete", "off_target"])
        summary = json.loads((root / "summary.json").read_text())
        self.assertEqual(summary["reads_total"], 2)
        self.assertEqual(summary["status_counts"]["complete"], 1)
        self.assertTrue((root / "products.tsv").is_file())

    def test_fastq_rejects_sequence_quality_length_mismatch(self):
        path = Path(self.tmp.name) / "bad.fastq"
        path.write_text("@bad\nACGT\n+\nIII\n", encoding="utf-8")
        with self.assertRaises(MoleculeAdmissionError):
            list(read_fastq(path))


if __name__ == "__main__":
    unittest.main()
