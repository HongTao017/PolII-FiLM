#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
py="$HOME/miniconda3/envs/groseq/bin/python"
plot_py=/mnt/HDD8TB/hongtao/Data_Seq/task/task_6/FIG5_SPLICE4_MLP_CHIP_ATTRIBUTION_20261013/analysis_venv/bin/python
exec 9>batch.lock
flock -n 9
printf '%s Waiting for validated FOXA2 pilot\n' "$(date -Is)"
flock pilot.lock true
"$py" -B - <<'PY'
import json
from pathlib import Path

expected = json.loads(Path("pilot_expected_usable.json").read_text())
pilot = json.loads(Path("pilot_status.json").read_text())
if pilot["status"] != "complete" or pilot["passed"] != len(expected):
    raise ValueError("FOXA2 pilot failed")
for accession, expected_count in expected.items():
    actual = json.loads((Path("counts") / f"{accession}.json").read_text())["usable"]
    if actual != expected_count:
        raise ValueError(f"{accession}: counted {actual} fragments; ENCODE reports {expected_count}")
print("PASS FOXA2 count denominators match ENCODE")
PY
"$plot_py" -B plot.py --pilot
printf '%s Starting full CPM counts with four workers\n' "$(date -Is)"
"$py" -B -u cpm.py run --workers 4
"$plot_py" -B plot.py
printf '%s Saved ChIP_fragment_CPM_vs_IG.pdf and summary.json\n' "$(date -Is)"
