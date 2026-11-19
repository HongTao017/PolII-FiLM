#!/usr/bin/env python3
"""Rename merged mNET-seq samples and synchronize sample IDs in PI outputs."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path


DEFAULT_PROJECT = Path("/mnt/HDD8TB/shared/Data_Seq/DATA/mnetseq_project")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=DEFAULT_PROJECT)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Apply changes. Without this option, only print a dry-run plan.",
    )
    return parser.parse_args()


def read_mapping(path: Path) -> tuple[dict[str, str], dict[str, dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"old_sample_id", "formal_sample_name", "gse", "gsm"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Mapping must contain: {', '.join(sorted(required))}")
    mapping = {row["old_sample_id"]: row["formal_sample_name"] for row in rows}
    metadata = {row["old_sample_id"]: row for row in rows}
    if len(mapping) != len(rows):
        raise ValueError("Duplicate old_sample_id values were found")
    if len(set(mapping.values())) != len(mapping):
        raise ValueError("Duplicate formal_sample_name values were found")
    if len(mapping) != 20:
        raise ValueError(f"Expected 20 mNET-seq samples, found {len(mapping)}")
    return mapping, metadata


def discover_one(project: Path, names: list[str], required: bool = True) -> Path | None:
    hits: list[Path] = []
    for name in names:
        hits.extend(path for path in project.rglob(name) if path.is_file())
    unique = sorted(set(hits))
    if len(unique) == 1:
        return unique[0]
    if not unique and not required:
        return None
    if not unique:
        raise FileNotFoundError(f"Could not find any of: {', '.join(names)}")
    raise ValueError(f"Multiple candidate files found: {', '.join(map(str, unique))}")


def prefer_or_discover(
    project: Path,
    preferred: list[Path],
    names: list[str],
    required: bool = True,
) -> Path | None:
    """Use a canonical project path before considering recursive search hits."""
    existing = [path for path in preferred if path.is_file()]
    if len(existing) == 1:
        return existing[0]
    if len(existing) > 1:
        raise ValueError(
            "Multiple canonical candidate files found: "
            + ", ".join(map(str, existing))
        )
    return discover_one(project, names, required=required)


def find_sample_bed_dir(project: Path) -> Path:
    candidates = [
        project / "results/mNET_3prime/sample_bed",
        project / "results/mNET_3prime/sample_merged",
        project / "rna_3prime_end/sample_merged",
    ]
    for path in candidates:
        if path.is_dir() and any(path.glob("*.bed.gz")):
            return path
    raise FileNotFoundError("Could not locate the merged mNET-seq BED directory")


def find_pi_dir(project: Path) -> Path:
    candidates = [project / "results/PI/mNET", project / "results/PI/mnet"]
    for path in candidates:
        if path.is_dir() and any(path.glob("*.mNET.PI.tsv")):
            return path
    raise FileNotFoundError("Could not locate the mNET-seq PI directory")


def planned_renames(directory: Path, mapping: dict[str, str]) -> list[tuple[Path, Path]]:
    moves: list[tuple[Path, Path]] = []
    for source in directory.iterdir():
        if not source.is_file():
            continue
        for old, new in mapping.items():
            if source.name == old or source.name.startswith(old + "."):
                target = source.with_name(new + source.name[len(old) :])
                moves.append((source, target))
                break
    return moves


def validate_moves(moves: list[tuple[Path, Path]]) -> None:
    targets = [target for _, target in moves]
    if len(targets) != len(set(targets)):
        raise ValueError("Two source files would be renamed to the same target")
    sources = {source for source, _ in moves}
    collisions = [target for _, target in moves if target.exists() and target not in sources]
    if collisions:
        raise FileExistsError("Rename targets already exist: " + ", ".join(map(str, collisions)))


def replace_path_sample(value: str, mapping: dict[str, str]) -> str:
    for old, new in mapping.items():
        if Path(value).name == old or Path(value).name.startswith(old + "."):
            return str(Path(value).with_name(new + Path(value).name[len(old) :]))
    return value


def rewrite_tsv(
    source: Path,
    output: Path,
    mapping: dict[str, str],
    metadata: dict[str, dict[str, str]],
) -> None:
    with source.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)

    if "sample_id" in fields and "original_sample_id" not in fields:
        fields.insert(fields.index("sample_id") + 1, "original_sample_id")
    if "sample_id" in fields and "formal_sample_name" not in fields:
        fields.insert(fields.index("sample_id") + 1, "formal_sample_name")

    path_columns = {
        "output_bed_gz",
        "output_3prime_bed",
        "output_pi",
    }
    for row in rows:
        old = row.get("sample_id", "")
        if old in mapping:
            row["original_sample_id"] = old
            row["formal_sample_name"] = mapping[old]
            row["sample_id"] = mapping[old]
            if "gse" in row:
                row["gse"] = metadata[old]["gse"]
            if "gsm" in row:
                row["gsm"] = metadata[old]["gsm"]
        for column in path_columns & row.keys():
            row[column] = replace_path_sample(row[column], mapping)

    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def rewrite_matrix(source: Path, output: Path, mapping: dict[str, str]) -> None:
    with source.open(encoding="utf-8-sig", newline="") as in_handle, output.open(
        "w", encoding="utf-8", newline=""
    ) as out_handle:
        reader = csv.reader(in_handle, delimiter="\t")
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
        header = next(reader)
        writer.writerow([header[0], *[mapping.get(value, value) for value in header[1:]]])
        writer.writerows(reader)


def main() -> int:
    args = parse_args()
    project = args.project.resolve()
    mapping, metadata = read_mapping(args.mapping.resolve())
    sample_bed_dir = find_sample_bed_dir(project)
    pi_dir = find_pi_dir(project)
    merge_status = prefer_or_discover(
        project,
        [
            project / "results/mnetseq_3prime_sample_merge_status.tsv",
            project / "results/mnetseq_sample_merge_status.tsv",
            project / "results/mNET_3prime/sample_merge_status.tsv",
        ],
        [
            "mnetseq_3prime_sample_merge_status.tsv",
            "mnetseq_sample_merge_status.tsv",
            "sample_merge_status.tsv",
        ],
    )
    pi_summary = prefer_or_discover(
        project,
        [project / "results/PI/mNET/pi_sample_summary.tsv"],
        ["pi_sample_summary.tsv"],
    )
    pi_matrix = prefer_or_discover(
        project,
        [project / "results/PI/mNET/mnetseq_pi_matrix.tsv"],
        ["mnetseq_pi_matrix.tsv"],
    )
    sample_audit = prefer_or_discover(
        project,
        [
            project / "results/mnetseq_sample_audit.tsv",
            project / "results/mnetseq_sample_final_audit.tsv",
        ],
        ["mnetseq_sample_audit.tsv", "mnetseq_sample_final_audit.tsv"],
        required=False,
    )

    moves = planned_renames(sample_bed_dir, mapping) + planned_renames(pi_dir, mapping)
    validate_moves(moves)
    found_old = {
        old
        for old in mapping
        if any(source.name == old or source.name.startswith(old + ".") for source, _ in moves)
    }
    missing = sorted(set(mapping) - found_old)
    if missing:
        raise FileNotFoundError("No merged or PI files found for: " + ", ".join(missing))

    tables = [merge_status, pi_summary, pi_matrix]
    if sample_audit:
        tables.append(sample_audit)
    print(f"Project: {project}")
    print(f"Merged BED directory: {sample_bed_dir}")
    print(f"PI directory: {pi_dir}")
    print(f"File renames: {len(moves)}")
    for source, target in moves:
        print(f"  {source.name} -> {target.name}")
    print("Tables to synchronize:")
    for table in tables:
        print(f"  {table}")

    if not args.apply:
        print("DRY RUN ONLY. Re-run with --apply after checking the plan.")
        return 0

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_dir = project / "results" / f"mnetseq_name_backup_{stamp}"
    backup_dir.mkdir(parents=True, exist_ok=False)
    manifest = {
        "created": stamp,
        "mapping": mapping,
        "moves": [[str(source), str(target)] for source, target in moves],
        "table_backups": {},
    }

    for table in tables:
        backup = backup_dir / table.name
        shutil.copy2(table, backup)
        manifest["table_backups"][str(table)] = str(backup)

    temporary_tables: list[tuple[Path, Path]] = []
    for table in tables:
        temp = table.with_name(table.name + f".rename_tmp_{stamp}")
        if table == pi_matrix:
            rewrite_matrix(table, temp, mapping)
        else:
            rewrite_tsv(table, temp, mapping, metadata)
        temporary_tables.append((temp, table))

    completed_moves: list[tuple[Path, Path]] = []
    try:
        for source, target in moves:
            os.replace(source, target)
            completed_moves.append((source, target))
        for temp, table in temporary_tables:
            os.replace(temp, table)
    except Exception:
        for temp, _ in temporary_tables:
            temp.unlink(missing_ok=True)
        for source, target in reversed(completed_moves):
            if target.exists() and not source.exists():
                os.replace(target, source)
        for original, backup in manifest["table_backups"].items():
            shutil.copy2(backup, original)
        raise

    rollback_manifest = backup_dir / "rollback_manifest.json"
    rollback_manifest.write_text(
        json.dumps(manifest, ensure_ascii=True, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Applied formal names to {len(mapping)} samples.")
    print(f"Backup and rollback manifest: {backup_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise SystemExit(1)
