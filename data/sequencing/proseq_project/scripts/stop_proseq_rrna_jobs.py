#!/usr/bin/env python3
"""Stop only this project's rRNA filtering processes."""

import os
import signal
import time
from pathlib import Path


PROJECT = "/mnt/HDD8TB/shared/Data_Seq/DATA/proseq_project"
TOKENS = ("filter_proseq_rrna.py", "human_rRNA", "rrna_filtered/final")


def matching_pids() -> list[int]:
    found = []
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        if pid == os.getpid():
            continue
        try:
            command = (entry / "cmdline").read_bytes().replace(b"\0", b" ").decode(
                "utf-8", "replace"
            )
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        if PROJECT in command and any(token in command for token in TOKENS):
            found.append(pid)
    return found


targets = matching_pids()
for pid in targets:
    try:
        os.kill(pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
print("SIGTERM:", " ".join(map(str, targets)) or "none", flush=True)

time.sleep(3)
remaining = matching_pids()
for pid in remaining:
    try:
        os.kill(pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
print("SIGKILL:", " ".join(map(str, remaining)) or "none", flush=True)
