#!/usr/bin/env python3
"""Create sanitized real-data fixtures from the lab's approved source archive."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path


SOURCE_SHA256 = "279ee3cb244436606e059bce49825b286b103b27c0af9355bd9852a6f0444029"
FASTQ_MEMBER_SHA256 = "e50a9b12241199b14be544bebe58442df94b1b2e76285f86aa47f92dd0d4a1e7"
BAM_MEMBER_SHA256 = "ebcb85a2d2c70f126987cfb77812407b27e8b50b056d2d934fbec6201d3e114b"
EXPECTED_FASTQ_RECORDS = 234
EXPECTED_BAM_RECORDS = 1329
EXPECTED_MAPPED_RECORDS = 1115
EXPECTED_REFERENCES = 48
REDISTRIBUTION_AUTHORIZATION_DATE = "2026-08-05"

FASTQ_SUFFIX = ".fastq.gz"
BAM_SUFFIX = ".bam"
ALIGNMENT_TAGS = {
    "AS", "MD", "NM", "SA", "cm", "de", "ms", "nn", "rl", "s1", "s2", "tp"
}


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def run(command: list[str], **kwargs) -> subprocess.CompletedProcess:
    kwargs.setdefault("encoding", "utf-8")
    kwargs.setdefault("errors", "replace")
    return subprocess.run(command, check=True, text=True, **kwargs)


def find_member(archive: zipfile.ZipFile, suffix: str) -> str:
    matches = [name for name in archive.namelist() if name.lower().endswith(suffix)]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one {suffix} member, found {len(matches)}")
    return matches[0]


def sanitize_fastq(source_data: bytes, destination: Path) -> dict[str, int]:
    record_count = 0
    base_count = 0
    source = gzip.GzipFile(fileobj=io.BytesIO(source_data), mode="rb")
    destination.parent.mkdir(parents=True, exist_ok=True)

    with source, destination.open("wb") as raw_output:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_output, mtime=0) as output:
            while True:
                lines = [source.readline() for _ in range(4)]
                if lines == [b"", b"", b"", b""]:
                    break
                if any(not line for line in lines):
                    raise RuntimeError("FASTQ ended inside a record")

                header, sequence, separator, quality = [line.rstrip(b"\r\n") for line in lines]
                if not header.startswith(b"@") or not separator.startswith(b"+"):
                    raise RuntimeError(f"Malformed FASTQ record {record_count + 1}")
                if len(sequence) != len(quality):
                    raise RuntimeError(f"Sequence/quality mismatch in record {record_count + 1}")

                record_count += 1
                base_count += len(sequence)
                read_id = f"fastq_test_read_{record_count:06d}".encode("ascii")
                sanitized_header = (
                    b"@" + read_id
                    + b"\tRG:Z:test_run_barcode19"
                    + b"\tSM:Z:barcode19"
                    + b"\tal:Z:barcode19"
                    + b"\tLB:Z:test_fixture"
                    + b"\tPU:Z:test_flowcell"
                )
                output.write(sanitized_header + b"\n")
                output.write(sequence + b"\n+\n" + quality + b"\n")

    if record_count != EXPECTED_FASTQ_RECORDS:
        raise RuntimeError(
            f"Expected {EXPECTED_FASTQ_RECORDS} FASTQ records, found {record_count}"
        )
    return {"records": record_count, "bases": base_count}


def sanitized_header_line(line: str) -> str | None:
    if line.startswith("@HD"):
        return "@HD\tVN:1.6\tSO:coordinate"
    if line.startswith("@SQ"):
        fields = line.rstrip("\n").split("\t")
        keep = [field for field in fields[1:] if field.startswith(("SN:", "LN:", "M5:"))]
        return "\t".join(["@SQ", *keep])
    if line.startswith(("@PG", "@RG", "@CO")):
        return None
    return line.rstrip("\n")


def sanitize_sam(source_bam: Path, destination_sam: Path, samtools: str) -> dict[str, int]:
    process = subprocess.Popen(
        [samtools, "view", "-h", str(source_bam)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
    )
    if process.stdout is None or process.stderr is None:
        raise RuntimeError("Could not capture samtools output")

    read_names: dict[str, str] = {}
    records = 0
    mapped = 0
    references = 0
    mapq_min: int | None = None
    mapq_max: int | None = None

    with destination_sam.open("w", encoding="utf-8", newline="\n") as output:
        inserted_read_group = False
        for line in process.stdout:
            if line.startswith("@"):
                header = sanitized_header_line(line)
                if header is not None:
                    output.write(header + "\n")
                    if header.startswith("@SQ"):
                        references += 1
                continue

            if not inserted_read_group:
                output.write(
                    "@RG\tID:test_run_barcode07\tSM:barcode07\tLB:test_fixture"
                    "\tPU:test_flowcell\tPL:ONT\n"
                )
                inserted_read_group = True

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise RuntimeError(f"Malformed SAM record {records + 1}")

            original_name = fields[0]
            if original_name not in read_names:
                read_names[original_name] = f"bam_test_read_{len(read_names) + 1:06d}"
            fields[0] = read_names[original_name]

            kept_tags = [tag for tag in fields[11:] if tag.split(":", 1)[0] in ALIGNMENT_TAGS]
            fields[11:] = ["RG:Z:test_run_barcode07", *kept_tags]
            output.write("\t".join(fields) + "\n")

            records += 1
            flag = int(fields[1])
            mapq = int(fields[4])
            if not flag & 4:
                mapped += 1
            mapq_min = mapq if mapq_min is None else min(mapq_min, mapq)
            mapq_max = mapq if mapq_max is None else max(mapq_max, mapq)

    stderr = process.stderr.read()
    returncode = process.wait()
    if returncode:
        raise RuntimeError(f"samtools view failed: {stderr.strip()}")
    if not inserted_read_group:
        raise RuntimeError("Source BAM contained no records")

    observed = {
        "records": records,
        "mapped_records": mapped,
        "unique_read_ids": len(read_names),
        "references": references,
        "mapq_min": mapq_min,
        "mapq_max": mapq_max,
    }
    expected = {
        "records": EXPECTED_BAM_RECORDS,
        "mapped_records": EXPECTED_MAPPED_RECORDS,
        "references": EXPECTED_REFERENCES,
    }
    for key, value in expected.items():
        if observed[key] != value:
            raise RuntimeError(f"Expected BAM {key}={value}, found {observed[key]}")
    return observed


def sanitize_bam(source_data: bytes, destination: Path, samtools: str) -> dict[str, int]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="fixture-bam-") as temporary:
        temporary_path = Path(temporary)
        source_bam = temporary_path / "source.bam"
        sanitized_sam = temporary_path / "sanitized.sam"
        source_bam.write_bytes(source_data)

        statistics = sanitize_sam(source_bam, sanitized_sam, samtools)
        run([
            samtools,
            "sort",
            "--no-PG",
            "-O",
            "BAM",
            "-o",
            str(destination),
            str(sanitized_sam),
        ])
        run([samtools, "index", str(destination)])
        run([samtools, "quickcheck", "-v", str(destination)])

        count = int(run([samtools, "view", "-c", str(destination)], capture_output=True).stdout)
        mapped = int(
            run([samtools, "view", "-c", "-F", "4", str(destination)], capture_output=True).stdout
        )
        if count != EXPECTED_BAM_RECORDS or mapped != EXPECTED_MAPPED_RECORDS:
            raise RuntimeError("Generated BAM record counts do not match the source fixture")

    return statistics


def write_manifest(
    destination: Path,
    fastq_path: Path,
    bam_path: Path,
    fastq_statistics: dict[str, int],
    bam_statistics: dict[str, int],
    samtools: str,
) -> None:
    version = run([samtools, "--version"], capture_output=True).stdout.splitlines()[0]
    manifest = {
        "schema_version": 1,
        "source": {
            "archive_name": "test_case.zip",
            "sha256": SOURCE_SHA256,
            "committed": False,
        },
        "provenance": {
            "data_type": "sanitized real human nanopore sequencing data",
            "fastq_and_bam_are_independent_runs": True,
            "sample_identifiers_replaced_with_anonymous_test_labels": True,
            "public_redistribution_authorized_by_lab": True,
            "redistribution_authorization_recorded_on": REDISTRIBUTION_AUTHORIZATION_DATE,
        },
        "sanitization": {
            "read_ids_replaced": True,
            "run_and_device_metadata_removed": True,
            "sequences_and_qualities_preserved": True,
            "bam_alignment_fields_preserved": True,
            "bam_tags_preserved": sorted(ALIGNMENT_TAGS | {"RG"}),
        },
        "files": {
            fastq_path.relative_to(destination).as_posix(): {
                "sha256": sha256_path(fastq_path),
                **fastq_statistics,
                "stage": "NanoTel and NanoTel filtration",
            },
            bam_path.relative_to(destination).as_posix(): {
                "sha256": sha256_path(bam_path),
                **bam_statistics,
                "stage": "aligned BAM summary and mapping",
            },
            bam_path.with_suffix(bam_path.suffix + ".bai").relative_to(destination).as_posix(): {
                "sha256": sha256_path(bam_path.with_suffix(bam_path.suffix + ".bai")),
            },
        },
        "preparation": {
            "script": "tests/fixtures/prepare_integration_fixtures.py",
            "samtools": version,
        },
    }
    (destination / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def prepare(source_archive: Path, destination: Path, samtools: str) -> None:
    if sha256_path(source_archive) != SOURCE_SHA256:
        raise RuntimeError(
            "Source archive checksum does not match the approved test_case.zip; refusing to continue"
        )
    if shutil.which(samtools) is None:
        raise RuntimeError(f"samtools executable not found: {samtools}")

    with zipfile.ZipFile(source_archive) as archive:
        fastq_member = find_member(archive, FASTQ_SUFFIX)
        bam_member = find_member(archive, BAM_SUFFIX)
        fastq_data = archive.read(fastq_member)
        bam_data = archive.read(bam_member)

    if sha256_bytes(fastq_data) != FASTQ_MEMBER_SHA256:
        raise RuntimeError("FASTQ member checksum does not match the approved source")
    if sha256_bytes(bam_data) != BAM_MEMBER_SHA256:
        raise RuntimeError("BAM member checksum does not match the approved source")

    fastq_path = destination / "fastq" / "barcode19" / "test_reads_barcode19.fastq.gz"
    bam_path = destination / "bam" / "barcode07" / "test_aligned_barcode07.bam"
    fastq_statistics = sanitize_fastq(fastq_data, fastq_path)
    bam_statistics = sanitize_bam(bam_data, bam_path, samtools)
    write_manifest(
        destination,
        fastq_path,
        bam_path,
        fastq_statistics,
        bam_statistics,
        samtools,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_archive", type=Path, help="Path to the ignored test_case.zip")
    parser.add_argument(
        "--destination",
        type=Path,
        default=Path(__file__).resolve().parent / "integration",
        help="Directory in which to create the sanitized fixtures",
    )
    parser.add_argument("--samtools", default="samtools", help="samtools executable")
    args = parser.parse_args()

    prepare(args.source_archive.resolve(), args.destination.resolve(), args.samtools)
    print(f"Sanitized integration fixtures written to {args.destination.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
