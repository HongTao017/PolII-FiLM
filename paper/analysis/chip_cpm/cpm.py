"""Replace the Fig5d ChIP horizontal coordinate with matched regional CPM."""
import argparse
import concurrent.futures as cf
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict

ROOT = Path(__file__).resolve().parent
DATA = Path("/mnt/HDD8TB/hongtao/Data_Seq")
ORIGINAL = DATA / "task/task_6/A_CHIP_ZEROCOUNT_GE40_1550_20261104"
SOURCE = DATA / "DATA/encode_seq/DNA/ChIP_seq/metadata/ENCODE_ChIP_seq_selected_manifest.tsv"
POINTS = ORIGINAL / "Fig5d_factor_1550genes/factor_source_cell_means.tsv"
HOST = "https://www.encodeproject.org"
ACCESSION = re.compile(r"^/files/(ENCFF[A-Z0-9]+)/$")
METADATA_FIELDS = ["accession", "derived_from"]
BAM_FIELDS = ["accession", "dataset", "assembly", "output_type", "status", "file_format",
              "file_size", "md5sum", "s3_uri", "href", "biological_replicates"]


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_json(path, data):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def log(message):
    print(message, flush=True)


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def request(url):
    return urllib.request.urlopen(
        urllib.request.Request(url, headers={"Accept": "application/json",
                                              "User-Agent": "PolII-ChIP-CPM/1.0"}),
        timeout=45,
    )


def search(filters, fields):
    query = [("type", "File"), ("format", "json"), ("limit", "all"), *filters]
    query.extend(("field", field) for field in fields)
    url = HOST + "/search/?" + urllib.parse.urlencode(query)
    for retry in range(4):
        try:
            with request(url) as response:
                result = json.load(response)
            found = result["@graph"]
            if len(found) != result["total"]:
                raise ValueError("Incomplete ENCODE search response")
            return found
        except urllib.error.HTTPError as error:
            if error.code == 404:
                return []
            if retry == 3 or (error.code < 500 and error.code != 429):
                raise
        except (OSError, ValueError):
            if retry == 3:
                raise
        time.sleep(2 ** retry)


def batch_search(key, keys, filters, fields):
    keys = sorted(set(keys))
    batches = [keys[i:i + 24] for i in range(0, len(keys), 24)]
    with cf.ThreadPoolExecutor(max_workers=5) as pool:
        futures = [pool.submit(search, filters + [(key, value) for value in batch], fields)
                   for batch in batches]
        return [record for future in cf.as_completed(futures) for record in future.result()]


def accession(path):
    if isinstance(path, dict):
        path = path.get("@id", "")
    match = ACCESSION.fullmatch(path)
    return match[1] if match else None


def geometry():
    genes = rows(ORIGINAL / "selected_genes.tsv")
    if len(genes) != 1550 or len({gene["gene"] for gene in genes}) != 1550:
        raise ValueError("Wrong original gene coverage")
    records = []
    for gene in genes:
        tss, tes = int(gene["tss_bed0"]), int(gene["tes_bed0"])
        if gene["strand"] == "+":
            pause = max(0, tss - 50), tss + 250
            body = tss + 250, tes + 1
        elif gene["strand"] == "-":
            pause = max(0, tss - 250), tss + 50
            body = tes, tss - 250 + 1
        else:
            raise ValueError("Unexpected strand")
        for label, (start, end) in (("pause", pause), ("gene_body", body)):
            records.append((gene["chrom"], start, max(start, end), gene["gene"] + "|" + label))
    return records


def write_bed():
    records = geometry()
    path = ROOT / "regions.bed"
    content = "".join(f"{chrom}\t{start}\t{end}\t{name}\n"
                      for chrom, start, end, name in records if start < end)
    if path.exists() and path.read_text() != content:
        raise ValueError("Frozen counting regions changed")
    path.write_text(content)
    return path, records


def prepare(pilot_only=False):
    original_points = rows(POINTS)
    if len(original_points) != 838 or len({r["pair_id"] for r in original_points}) != 838:
        raise ValueError("Wrong original group coverage")

    def key(cell, target):
        return tuple(re.sub(r"[^a-z0-9]", "", s.lower()) for s in (cell, target))

    selected = defaultdict(list)
    for entry in rows(SOURCE):
        selected[key(entry["cell_line"], entry["target"])].append(entry)

    groups = []
    for point in original_points:
        if pilot_only and point["pair_id"] != "HepG2_FOXA2":
            continue
        (peak,) = selected[key(point["cell"], point["tf"])]
        groups.append({**point, "experiment": peak["experiment_accession"],
                       "peak": peak["file_accession"],
                       "replicates": sorted(int(v) for v in peak["biological_replicates"].split(",") if v)})
    by_dataset = defaultdict(dict)
    candidates = batch_search("dataset", [f'/experiments/{g["experiment"]}/' for g in groups],
                              [("file_format", "bam"), ("output_type", "alignments"),
                               ("status", "released"), ("assembly", "GRCh38")], BAM_FIELDS)
    for bam in candidates:
        by_dataset[bam["dataset"]][bam["accession"]] = bam
    log(f"ENCODE candidate filtered GRCh38 BAMs: {len(candidates)}")

    pending = [(group, {group["peak"]}, set(), defaultdict(set)) for group in groups]
    cache = {}
    for depth in range(7):
        needed = {item for _, front, _, _ in pending for item in front} - cache.keys()
        cache.update({record["accession"]: record.get("derived_from", [])
                      for record in batch_search("accession", needed, [], METADATA_FIELDS)})
        for item in needed:
            cache.setdefault(item, [])
        unresolved = []
        for group, front, seen, chosen in pending:
            bam_map = by_dataset[f'/experiments/{group["experiment"]}/']
            found = defaultdict(set)
            for item in front:
                if item in bam_map:
                    reps = bam_map[item].get("biological_replicates", [])
                    if len(reps) == 1 and reps[0] in group["replicates"]:
                        found[reps[0]].add(item)
            for rep, ids in found.items():
                if rep not in chosen:
                    chosen[rep] = ids
            next_front = {parent for item in front if item not in bam_map
                          for path in cache[item] if (parent := accession(path))}
            seen.update(front)
            next_front -= seen
            if next_front and set(chosen) != set(group["replicates"]) and depth < 6:
                unresolved.append((group, next_front, seen, chosen))
                continue
            group["bam_matches"] = {str(rep): sorted(ids) for rep, ids in sorted(chosen.items())}
            if group["replicates"] and set(chosen) == set(group["replicates"]) and all(len(ids) == 1 for ids in chosen.values()):
                group["bams"] = [bam_map[next(iter(chosen[rep]))] for rep in group["replicates"]]
                group["status"] = "matched"
            else:
                group["bams"] = []
                group["status"] = "incomplete_or_ambiguous_peak_ancestry"
        pending = unresolved
        log(f"Ancestry depth {depth}: {len(pending)} groups still searching")
        if not pending:
            break
    if len(groups) != (1 if pilot_only else 838):
        raise ValueError("Wrong planned group count")
    beds, regions = write_bed()
    out = ROOT / ("pilot_manifest.json" if pilot_only else "manifest.json")
    distinct = {bam["accession"]: bam for g in groups for bam in g["bams"]}
    stats = dict(groups=len(groups), matched_groups=sum(g["status"] == "matched" for g in groups),
                 bams=len(distinct), bytes=sum(int(b["file_size"]) for b in distinct.values()),
                 empty_regions=sum(start == end for _, start, end, _ in regions))
    atomic_json(out, dict(stats=stats, source_sha256=digest(SOURCE), bed_sha256=digest(beds),
                          points_sha256=digest(POINTS), groups=groups))
    log(stats)


def bam_url(bam, retry):
    s3 = bam.get("s3_uri", "")
    if retry != 1 and s3.startswith("s3://encode-public/"):
        return "https://encode-public.s3.amazonaws.com/" + s3.removeprefix("s3://encode-public/")
    return HOST + bam["href"]


def count_bam(bam, bed, region_names, out_dir):
    code = bam["accession"]
    result = out_dir / f"{code}.json"
    if result.exists():
        cached = json.loads(result.read_text())
        if (cached["md5sum"] == bam["md5sum"] and cached["bed_sha256"] == digest(bed)
                and cached["counter_sha256"] == digest(ROOT / "count_fragments")
                and set(cached["counts"]) == region_names):
            return code
        raise ValueError(f"{code}: saved count provenance mismatch")
    for retry in range(3):
        try:
            log(f"Counting {code}, {bam['file_size'] / 1e9:.2f} GB (attempt {retry + 1})")
            with tempfile.TemporaryFile() as stdout, tempfile.TemporaryFile() as stderr:
                cmd = [str(ROOT / "count_fragments"), str(bed)]
                with subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=stdout, stderr=stderr) as counter:
                    count_bytes = 0
                    checksum = hashlib.md5()
                    try:
                        with urllib.request.urlopen(bam_url(bam, retry), timeout=60) as download:
                            while block := download.read(1024 * 1024):
                                checksum.update(block)
                                count_bytes += len(block)
                                counter.stdin.write(block)
                    finally:
                        counter.stdin.close()
                    status = counter.wait()
                stdout.seek(0)
                stderr.seek(0)
                error = stderr.read().decode(errors="replace").strip()
                if status or count_bytes != bam["file_size"] or checksum.hexdigest() != bam["md5sum"]:
                    raise ValueError(f"{code}: corrupt or rejected BAM ({status}, {count_bytes} bytes): {error}")
                counts = {}
                usable = 0
                for line in stdout.read().decode().splitlines():
                    label, count, *rest = line.split("\t")
                    if label == "#usable":
                        usable, layout = int(count), rest[0]
                    else:
                        counts[label] = int(count)
                if usable <= 0 or set(counts) != region_names or any(not 0 <= n <= usable for n in counts.values()):
                    raise ValueError(f"{code}: incorrect region counts")
                atomic_json(result, dict(accession=code, md5sum=checksum.hexdigest(), usable=usable,
                                         layout=layout, bed_sha256=digest(bed),
                                         counter_sha256=digest(ROOT / "count_fragments"), counts=counts))
            log(f"PASS {code}: {usable:,} usable {layout}")
            return code
        except (OSError, ValueError) as error:
            log(f"{code}: {error}")
            if retry == 2:
                raise
            time.sleep(5 * (retry + 1))


def run(workers, pilot_only=False):
    manifest = json.loads((ROOT / ("pilot_manifest.json" if pilot_only else "manifest.json")).read_text())
    bed, regions = write_bed()
    if digest(bed) != manifest["bed_sha256"]:
        raise ValueError("Changed region coordinates")
    bams = {bam["accession"]: bam for group in manifest["groups"] for bam in group["bams"]}
    if not bams:
        raise ValueError("No matched BAMs")
    out_dir = ROOT / "counts"
    out_dir.mkdir(exist_ok=True)
    region_names = {name for _, start, end, name in regions if start < end}
    failures = {}
    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(count_bam, bam, bed, region_names, out_dir): code for code, bam in bams.items()}
        for future in cf.as_completed(futures):
            try:
                future.result()
            except (OSError, ValueError) as error:
                failures[futures[future]] = str(error)
                log(f"FAILED {futures[future]}: {error}")
    atomic_json(ROOT / ("pilot_status.json" if pilot_only else "status.json"),
                dict(status="complete" if not failures else "counting_failed", total=len(bams),
                     passed=len(bams) - len(failures), failed=failures))
    if failures:
        raise RuntimeError(f"{len(failures)} BAMs failed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "run"))
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--workers", type=int, choices=range(1, 7), default=3)
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.pilot)
    else:
        run(args.workers, args.pilot)


if __name__ == "__main__":
    main()
