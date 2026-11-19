#!/usr/bin/env python3
"""Split a mixed download manifest into true mNET-seq and excluded assay rows."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


NON_MNET_RUNS = {"SRR14777732", "SRR14777733", "SRR14777736", "SRR14777737"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--mnet-output", required=True, type=Path)
    parser.add_argument("--excluded-output", required=True, type=Path)
    args = parser.parse_args()

    with args.manifest.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)

    mnet = [row for row in rows if row["run"] not in NON_MNET_RUNS]
    excluded = [row for row in rows if row["run"] in NON_MNET_RUNS]
    for path, selected in ((args.mnet_output, mnet), (args.excluded_output, excluded)):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(selected)

    print(
        f"mNET_rows={len(mnet)} mNET_runs={len({row['run'] for row in mnet})} "
        f"excluded_rows={len(excluded)} excluded_runs={len({row['run'] for row in excluded})}"
    )


if __name__ == "__main__":
    main()
