import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import dorado_workflow
from dorado_workflow.processors.aligner import AlignmentProcessor
from dorado_workflow.processors.r_analyzer import RAnalyzer


RUN_REAL_INTEGRATION_TESTS = os.environ.get("RUN_REAL_INTEGRATION_TESTS") == "1"
REAL_INTEGRATION_OUTPUT_DIR = os.environ.get("REAL_INTEGRATION_OUTPUT_DIR")
FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures" / "integration"
EXPECTED_RESULTS_PATH = FIXTURE_ROOT / "expected_results.json"
FASTQ_PATH = FIXTURE_ROOT / "fastq" / "barcode19" / "test_reads_barcode19.fastq.gz"
BAM_PATH = FIXTURE_ROOT / "bam" / "barcode07" / "test_aligned_barcode07.bam"
EXPECTED_SUMMARY_COLUMNS = {
    "sequence_ID",
    "sequence_length",
    "telo_density_mismatch",
    "Telomere_start_mismatch",
    "Telomere_length_mismatch",
}


def run(command, **kwargs):
    completed = subprocess.run(
        [str(part) for part in command],
        check=False,
        text=True,
        capture_output=True,
        **kwargs,
    )
    if completed.returncode:
        rendered = " ".join(str(part) for part in command)
        raise AssertionError(
            f"Command failed with exit code {completed.returncode}: {rendered}\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    return completed


def require_tool(name):
    path = shutil.which(name)
    if path is None:
        raise AssertionError(f"{name} is required when RUN_REAL_INTEGRATION_TESTS=1")
    return path


def read_csv_rows(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


@unittest.skipUnless(
    RUN_REAL_INTEGRATION_TESTS,
    "set RUN_REAL_INTEGRATION_TESTS=1 to run sanitized real-data integration tests",
)
class RealDataIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.samtools = require_tool("samtools")
        cls.rscript = require_tool("Rscript")
        cls.workflow_root = Path(dorado_workflow.__file__).resolve().parent
        cls.baseline = json.loads(EXPECTED_RESULTS_PATH.read_text(encoding="utf-8"))
        if REAL_INTEGRATION_OUTPUT_DIR:
            cls.temporary = None
            cls.work = Path(REAL_INTEGRATION_OUTPUT_DIR).expanduser().resolve()
            cls.work.mkdir(parents=True, exist_ok=False)
            print(f"Persistent integration outputs: {cls.work}")
        else:
            cls.temporary = tempfile.TemporaryDirectory(prefix="telomere-real-integration-")
            cls.work = Path(cls.temporary.name)

    @classmethod
    def tearDownClass(cls):
        if cls.temporary is not None:
            cls.temporary.cleanup()

    def test_fastq_integrity_and_sanitization(self):
        self.assertTrue(FASTQ_PATH.is_file())
        records = 0
        with gzip.open(FASTQ_PATH, "rt", encoding="ascii") as handle:
            while True:
                lines = [handle.readline().rstrip("\r\n") for _ in range(4)]
                if lines == ["", "", "", ""]:
                    break
                records += 1
                self.assertRegex(lines[0], rf"^@fastq_test_read_{records:06d}\t")
                self.assertIn("RG:Z:test_run_barcode19", lines[0])
                self.assertIn("SM:Z:barcode19", lines[0])
                self.assertIn("LB:Z:test_fixture", lines[0])
                self.assertIn("PU:Z:test_flowcell", lines[0])
                self.assertEqual(lines[2], "+")
                self.assertEqual(len(lines[1]), len(lines[3]))
                self.assertNotRegex(lines[0], r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}")
        self.assertEqual(records, self.baseline["expected"]["fastq"]["input_records"])

    def test_fastq_nanotel_and_filtration(self):
        nanotel_output = self.work / "nanotel_raw" / "barcode19"
        nanotel_output.mkdir(parents=True)
        nanotel_script = self.workflow_root / "data" / "NanoTel.R"
        nanotel_parameters = self.baseline["parameters"]["nanotel"]
        run([
            self.rscript,
            "--vanilla",
            nanotel_script,
            "--input_path",
            FASTQ_PATH.parent,
            "--save_path",
            nanotel_output,
            "--summary_only",
            "--analysis",
            "--patterns",
            nanotel_parameters["patterns"],
        ])

        summary_file = nanotel_output / "barcode19_summary.csv"
        read_id_files = list(nanotel_output.glob("*_reads_ids.txt"))
        log_files = list(nanotel_output.rglob("*.log"))
        self.assertTrue(summary_file.is_file(), "NanoTel did not create its raw summary CSV")
        self.assertTrue(read_id_files, "NanoTel did not create a read-ID file")
        self.assertTrue(log_files, "NanoTel did not create a log")

        rows = read_csv_rows(summary_file)
        self.assertTrue(rows, "NanoTel summary is empty")
        self.assertTrue(EXPECTED_SUMMARY_COLUMNS.issubset(rows[0]))
        observed_ids = {row["sequence_ID"].split()[0] for row in rows}
        expected_fastq = self.baseline["expected"]["fastq"]
        self.assertEqual(observed_ids, set(expected_fastq["detected_read_ids"]))
        for row in rows:
            self.assertRegex(row["sequence_ID"].split()[0], r"^fastq_test_read_\d{6}$")
            for column in EXPECTED_SUMMARY_COLUMNS - {"sequence_ID"}:
                self.assertGreaterEqual(float(row[column]), 0)

        analysis_expected = expected_fastq["analysis"]
        analysis_summary = nanotel_output / "barcode19_filtered_sorted_summary.csv"
        analysis_results = nanotel_output / "barcode19_results.txt"
        analysis_plot = nanotel_output / "barcode19_telomere_plot.png"

        self.assertTrue(analysis_summary.is_file(), "NanoTel analysis summary is missing")
        analysis_rows = read_csv_rows(analysis_summary)
        self.assertEqual(len(analysis_rows), analysis_expected["filtered_records"])
        self.assertEqual(
            analysis_rows[0]["sequence_ID"].split()[0],
            analysis_expected["filtered_read_id"],
        )
        for column in ("TelLenMM_RunningMed", "SeqLen_minus_RunMed"):
            self.assertIn(column, analysis_rows[0])
            self.assertTrue(float(analysis_rows[0][column]) >= 0)

        self.assertTrue(analysis_results.is_file(), "NanoTel analysis results are missing")
        results_text = analysis_results.read_text(encoding="utf-8")
        self.assertIn(
            "Number of telomeric reads (post-filtration) : "
            f"{analysis_expected['filtered_records']}",
            results_text,
        )
        self.assertIn("Median Telomeric Length (post-filtration)", results_text)
        if not analysis_expected["km_metrics_enabled"]:
            self.assertNotIn("KM Median", results_text)
            self.assertNotIn("Expected KM Median Bias", results_text)

        self.assertTrue(analysis_plot.is_file(), "NanoTel analysis plot is missing")
        self.assertGreater(analysis_plot.stat().st_size, 1024)
        self.assertEqual(analysis_plot.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")

        filtered_output = self.work / "nanotel_filtered"
        config_path = self.work / "nanotel_filter_config.json"
        filtration_parameters = self.baseline["parameters"]["filtration"]
        config_path.write_text(
            json.dumps({
                "input_dir": str(nanotel_output.parent),
                "output_dir": str(filtered_output),
                "density_threshold": filtration_parameters["density_threshold"],
                "max_telomere_start": filtration_parameters["max_telomere_start"],
            }),
            encoding="utf-8",
        )
        run([
            self.rscript,
            "--vanilla",
            self.workflow_root / "r_analysis" / "batch_nanotel_analysis.R",
            config_path,
        ])
        filtered_files = list(filtered_output.rglob("filtered_summary*.csv"))
        self.assertTrue(filtered_files, "NanoTel filtration produced no summary")
        filtered_rows = read_csv_rows(filtered_files[0])
        self.assertEqual(len(filtered_rows), expected_fastq["filtered_records"])
        self.assertIn("read_id", filtered_rows[0])
        filtered = filtered_rows[0]
        self.assertEqual(filtered["read_id"].split()[0], expected_fastq["filtered_read_id"])
        self.assertEqual(int(filtered["sequence_length"]), expected_fastq["filtered_sequence_length"])
        self.assertEqual(
            int(filtered["Telomere_length_mismatch"]),
            expected_fastq["filtered_telomere_length_mismatch"],
        )
        self.assertEqual(
            int(filtered["Telomere_start_mismatch"]),
            expected_fastq["filtered_telomere_start_mismatch"],
        )
        self.assertAlmostEqual(
            float(filtered["telo_density_mismatch"]),
            expected_fastq["filtered_telomere_density_mismatch"],
            places=12,
        )

    def test_bam_integrity_sanitization_mapping_and_no_methylation(self):
        self.assertTrue(BAM_PATH.is_file())
        self.assertTrue(Path(str(BAM_PATH) + ".bai").is_file())
        run([self.samtools, "quickcheck", "-v", BAM_PATH])

        header = run([self.samtools, "view", "--no-PG", "-H", BAM_PATH]).stdout
        expected_bam = self.baseline["expected"]["bam"]
        self.assertEqual(
            sum(line.startswith("@SQ") for line in header.splitlines()),
            expected_bam["references"],
        )
        self.assertIn("ID:test_run_barcode07", header)
        for forbidden in ("/home/", "tzfati", "PBK33040", "RTX 4090", "2026-05"):
            self.assertNotIn(forbidden, header)

        records = 0
        mapped = 0
        mapqs = []
        has_modified_base_tags = False
        view = subprocess.Popen(
            [self.samtools, "view", str(BAM_PATH)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        self.assertIsNotNone(view.stdout)
        for line in view.stdout:
            fields = line.rstrip("\n").split("\t")
            records += 1
            self.assertRegex(fields[0], r"^bam_test_read_\d{6}$")
            self.assertIn("RG:Z:test_run_barcode07", fields[11:])
            flag = int(fields[1])
            mapped += not bool(flag & 4)
            mapqs.append(int(fields[4]))
            tag_names = {field.split(":", 1)[0] for field in fields[11:]}
            has_modified_base_tags |= bool({"MM", "ML"} & tag_names)
        stderr = view.stderr.read() if view.stderr else ""
        self.assertEqual(view.wait(), 0, stderr)
        if view.stdout:
            view.stdout.close()
        if view.stderr:
            view.stderr.close()
        self.assertEqual(records, expected_bam["records"])
        self.assertEqual(mapped, expected_bam["mapped_records"])
        self.assertEqual(records - mapped, expected_bam["unmapped_records"])
        self.assertEqual((min(mapqs), max(mapqs)), (expected_bam["mapq_min"], expected_bam["mapq_max"]))
        self.assertEqual(has_modified_base_tags, expected_bam["has_modified_base_tags"])

        alignment_dir = self.work / "aligned"
        alignment_dir.mkdir(exist_ok=True)
        logger = SimpleNamespace(info=lambda *args, **kwargs: None)
        aligner = object.__new__(AlignmentProcessor)
        aligner.output_dir = alignment_dir
        aligner.context = SimpleNamespace(logger=logger)
        summary_path = aligner._write_alignment_summary([BAM_PATH])

        with summary_path.open(encoding="utf-8") as handle:
            summary_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(
            list(summary_rows[0]),
            [
                "read_id",
                "alignment_direction",
                "alignment_genome",
                "alignment_genome_start",
                "alignment_genome_end",
                "alignment_mapq",
            ],
        )
        expected_mapping = self.baseline["expected"]["mapping"]
        mapping_parameters = self.baseline["parameters"]["mapping"]
        self.assertEqual(len(summary_rows), expected_mapping["alignment_summary_rows"])
        self.assertEqual(
            sum(int(row["alignment_mapq"]) >= mapping_parameters["min_mapq"] for row in summary_rows),
            expected_bam["mapq_at_least_minimum"],
        )

        eligible = []
        for row in summary_rows:
            if row["alignment_direction"] == "*" or int(row["alignment_mapq"]) < mapping_parameters["min_mapq"]:
                continue
            reference = row["alignment_genome"].lower()
            start = int(row["alignment_genome_start"])
            end = int(row["alignment_genome_end"])
            if (
                "head" in reference and start <= mapping_parameters["head_max_start"]
            ) or (
                "tail" in reference and end >= mapping_parameters["tail_min_end"]
            ):
                eligible.append(row["read_id"])
        self.assertEqual(len(eligible), expected_mapping["after_position_filtering"])

        bam_dir = self.work / "bam" / "barcode07"
        bam_dir.mkdir(parents=True)
        staged_bam = bam_dir / BAM_PATH.name
        shutil.copy2(BAM_PATH, staged_bam)
        shutil.copy2(Path(str(BAM_PATH) + ".bai"), Path(str(staged_bam) + ".bai"))

        methylation_warnings = []
        analyzer = object.__new__(RAnalyzer)
        analyzer.context = SimpleNamespace(
            path_manager=SimpleNamespace(get_aligned_dir_path=lambda: bam_dir.parent),
            logger=SimpleNamespace(
                warning=methylation_warnings.append,
                info=lambda *args, **kwargs: None,
            ),
        )
        self.assertTrue(analyzer._check_bam_methylation(read_limit=1000))
        self.assertTrue(
            any("No methylation data found" in message for message in methylation_warnings)
        )

        filtered_dir = self.work / "mapping_input" / "barcode07"
        filtered_dir.mkdir(parents=True)
        mapping_input = filtered_dir / "filtered_summaryBARCODE07.csv"
        with mapping_input.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["read_id", "sequence_length"])
            writer.writeheader()
            for read_id in dict.fromkeys(eligible):
                writer.writerow({"read_id": read_id, "sequence_length": 0})

        mapping_output = self.work / "mapping_output"
        mapping_config = self.work / "mapping_config.json"
        mapping_config.write_text(
            json.dumps({
                "alignment_summary_path": str(summary_path),
                "filtered_nanotel_dir": str(filtered_dir.parent),
                "bam_dir": str(bam_dir.parent),
                "output_dir": str(mapping_output),
                "min_mapq": mapping_parameters["min_mapq"],
                "head_max_start": mapping_parameters["head_max_start"],
                "tail_min_end": mapping_parameters["tail_min_end"],
                "run_modkit_pileup": mapping_parameters["run_modkit_pileup"],
            }),
            encoding="utf-8",
        )
        run([
            self.rscript,
            "--vanilla",
            self.workflow_root / "r_analysis" / "batch_mapping_analysis.R",
            mapping_config,
        ])

        mapped_csvs = list(mapping_output.rglob("mapped*_combined.csv"))
        filtered_bams = list(mapping_output.rglob("filtered_*.bam"))
        self.assertTrue(mapped_csvs)
        self.assertEqual(len(read_csv_rows(mapped_csvs[0])), expected_mapping["combined_csv_rows"])
        self.assertTrue(filtered_bams)
        self.assertTrue(Path(str(filtered_bams[0]) + ".bai").is_file())
        run([self.samtools, "quickcheck", "-v", filtered_bams[0]])
        filtered_count = int(run([self.samtools, "view", "-c", filtered_bams[0]]).stdout)
        self.assertEqual(filtered_count, expected_mapping["filtered_bam_records"])
        self.assertEqual(
            len(list(mapping_output.rglob("*.bed"))),
            expected_mapping["methylation_bed_files"],
        )

        diagnostic_files = list(mapping_output.rglob("mapping_diagnostics_*.csv"))
        self.assertEqual(len(diagnostic_files), 1)
        observed_diagnostics = {
            row["step"]: int(row["count"])
            for row in read_csv_rows(diagnostic_files[0])
        }
        expected_diagnostics = {
            key: expected_mapping[key]
            for key in (
                "alignment_summary_rows",
                "nanotel_filtered_rows",
                "after_nanotel_read_id_match",
                "after_removing_unmapped",
                "after_mapq_filtering",
                "after_position_filtering",
                "after_strand_position_filtering",
            )
        }
        self.assertEqual(observed_diagnostics, expected_diagnostics)


if __name__ == "__main__":
    unittest.main()
