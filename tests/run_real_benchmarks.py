#!/usr/bin/env python3
"""
Real-genome assembly benchmark and SHA-256 bit-identity verification suite.
Compares frozen unmodified hifiasm_vanilla (reading staged disk files and writing *.bin checkpoints)
against single-pass streaming hifiasm across:
  1. Pure-RAM StreamStore (--stream-ram-gb 64 --no-write-bin, compressed FIFO stream)
  2. Pure-RAM StreamStore + parallel upstream producer (pigz -dc | hifiasm - --no-write-bin)
  3. Pure O_TMPFILE Spool (--stream-ram-gb 0 --no-write-bin, zero persistent StreamStore RAM)
  4. Hybrid RAM + O_TMPFILE Spool (--stream-ram-gb <half_store> --no-write-bin)
Records wall-clock runtime, Pass 1+2 I/O time, /usr/bin/time -v Peak RSS, disk bytes written,
GFA assembly contiguity (total bp, contig count, N50), and 100% SHA-256 checksum identity.
"""

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def gfa_contiguity_stats(gfa_path: Path) -> dict:
    lengths = []
    if not gfa_path.exists():
        return {"n_contigs": 0, "total_bp": 0, "n50_bp": 0, "max_bp": 0}
    with open(gfa_path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("S\t"):
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3 and parts[2] != "*":
                    lengths.append(len(parts[2]))
                else:
                    for tag in parts[3:]:
                        if tag.startswith("LN:i:"):
                            lengths.append(int(tag[5:]))
                            break
    if not lengths:
        return {"n_contigs": 0, "total_bp": 0, "n50_bp": 0, "max_bp": 0}
    lengths.sort(reverse=True)
    total = sum(lengths)
    half = total / 2.0
    acc = 0
    n50 = 0
    for l in lengths:
        acc += l
        if acc >= half:
            n50 = l
            break
    return {
        "n_contigs": len(lengths),
        "total_bp": total,
        "n50_bp": n50,
        "max_bp": lengths[0],
    }


def collect_outputs(prefix: Path) -> dict:
    parent = prefix.parent
    base = prefix.name
    outputs = {}
    for p in sorted(parent.glob(f"{base}.*")):
        if p.name.endswith(".gfa") or p.name.endswith(".bed"):
            rel = p.name[len(base) + 1 :]
            entry = {
                "sha256": sha256_file(p),
                "bytes": p.stat().st_size,
            }
            if p.name.endswith(".gfa") and not p.name.endswith(".noseq.gfa"):
                entry["contiguity"] = gfa_contiguity_stats(p)
            outputs[rel] = entry
    return outputs


def measure_disk_artifacts(prefix: Path) -> dict:
    parent = prefix.parent
    base = prefix.name
    bin_bytes = 0
    bin_files = []
    gfa_bed_bytes = 0
    for p in sorted(parent.glob(f"{base}.*")):
        sz = p.stat().st_size
        if p.name.endswith(".gfa") or p.name.endswith(".bed"):
            gfa_bed_bytes += sz
        else:
            bin_bytes += sz
            bin_files.append({"name": p.name[len(base) + 1 :], "bytes": sz})
    return {
        "checkpoint_bin_count": len(bin_files),
        "checkpoint_bin_bytes": bin_bytes,
        "checkpoint_bin_gb": round(bin_bytes / (1024 ** 3), 4),
        "gfa_bed_bytes": gfa_bed_bytes,
        "gfa_bed_gb": round(gfa_bed_bytes / (1024 ** 3), 4),
    }


def parse_log_metrics(log_text: str, time_text: str) -> dict:
    rss_kb = 0
    m_rss = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", time_text)
    if m_rss:
        rss_kb = int(m_rss.group(1))

    fs_out = 0
    m_fso = re.search(r"File system outputs:\s*(\d+)", time_text)
    if m_fso:
        fs_out = int(m_fso.group(1))

    # Extract hifiasm internal timestamps for Pass 1 (ha_ft_gen) and Pass 2 (first ha_pt_gen)
    ft_sec = None
    m_ft = re.search(r"\[M::ha_ft_gen\]\s*Real time:\s*([0-9.]+)\s*sec", log_text)
    if m_ft:
        ft_sec = float(m_ft.group(1))
    pt_matches = re.findall(r"\[M::ha_pt_gen\]\s*Real time:\s*([0-9.]+)\s*sec", log_text)
    pt0_sec = float(pt_matches[0]) if pt_matches else None

    st_info = None
    m_st = re.search(
        r"StreamStore captured (\d+) reads \((\d+) bp\) peak_ram=([0-9.]+) GB spilled_disk=([0-9.]+) GB",
        log_text,
    )
    if m_st:
        st_info = {
            "reads": int(m_st.group(1)),
            "bp": int(m_st.group(2)),
            "peak_ram_gb": float(m_st.group(3)),
            "spilled_disk_gb": float(m_st.group(4)),
        }

    return {
        "peak_rss_kb": rss_kb,
        "peak_rss_gb": round(rss_kb / 1048576.0, 4),
        "fs_outputs_blocks": fs_out,
        "pass1_ha_ft_sec": ft_sec,
        "pass2_ha_pt0_sec": pt0_sec,
        "pass12_sum_sec": round((ft_sec or 0.0) + (pt0_sec or 0.0), 3) if (ft_sec is not None and pt0_sec is not None) else None,
        "stream_store": st_info,
    }


def run_timed_vanilla(exe: Path, args: list, work_dir: Path, env: dict) -> tuple:
    work_dir.mkdir(parents=True, exist_ok=True)
    log_path = work_dir / "hifiasm.log"
    time_path = work_dir / "usr_time.txt"
    cmd = ["/usr/bin/time", "-v", "-o", str(time_path), str(exe)] + args
    t0 = time.perf_counter()
    with open(log_path, "wb") as flog:
        proc = subprocess.run(cmd, stdout=flog, stderr=subprocess.STDOUT, env=env)
    wall = time.perf_counter() - t0
    if proc.returncode != 0:
        tail = log_path.read_text(errors="replace")[-3000:]
        raise RuntimeError(f"Vanilla run failed (rc={proc.returncode}):\n{tail}")
    return wall, log_path.read_text(errors="replace"), time_path.read_text(errors="replace")


def run_timed_stream(
    exe: Path,
    args_before_reads: list,
    read_src: Path,
    work_dir: Path,
    env: dict,
    decompress_upstream: bool = False,
    hic_srcs: tuple = None,
) -> tuple:
    work_dir.mkdir(parents=True, exist_ok=True)
    log_path = work_dir / "hifiasm.log"
    time_path = work_dir / "usr_time.txt"
    writers = []
    cmd = ["/usr/bin/time", "-v", "-o", str(time_path), str(exe)] + list(args_before_reads)

    if hic_srcs is not None:
        fifo_h1 = work_dir / "hic_R1.fifo"
        fifo_h2 = work_dir / "hic_R2.fifo"
        os.mkfifo(fifo_h1)
        os.mkfifo(fifo_h2)
        if decompress_upstream:
            c1 = f'pigz -dc -p 3 "{hic_srcs[0]}" > "{fifo_h1}"'
            c2 = f'pigz -dc -p 3 "{hic_srcs[1]}" > "{fifo_h2}"'
        else:
            c1 = f'cat "{hic_srcs[0]}" > "{fifo_h1}"'
            c2 = f'cat "{hic_srcs[1]}" > "{fifo_h2}"'
        writers.append(subprocess.Popen(["bash", "-c", c1], env=env))
        writers.append(subprocess.Popen(["bash", "-c", c2], env=env))
        cmd.extend(["--h1", str(fifo_h1), "--h2", str(fifo_h2)])

    fifo_reads = work_dir / "long_reads.fifo"
    os.mkfifo(fifo_reads)
    if decompress_upstream:
        cr = f'pigz -dc -p 6 "{read_src}" > "{fifo_reads}"'
    else:
        cr = f'cat "{read_src}" > "{fifo_reads}"'
    writers.append(subprocess.Popen(["bash", "-c", cr], env=env))
    cmd.append("-")

    t0 = time.perf_counter()
    with open(fifo_reads, "rb") as fin, open(log_path, "wb") as flog:
        proc = subprocess.Popen(cmd, stdin=fin, stdout=flog, stderr=subprocess.STDOUT, env=env)
        rc = proc.wait()
    for w in writers:
        wrc = w.wait()
        if wrc != 0:
            raise RuntimeError(f"Upstream stream producer failed with rc={wrc}")
    wall = time.perf_counter() - t0
    if rc != 0:
        tail = log_path.read_text(errors="replace")[-3000:]
        raise RuntimeError(f"Stream run failed (rc={rc}):\n{tail}")
    return wall, log_path.read_text(errors="replace"), time_path.read_text(errors="replace")


def run_dataset_benchmark(
    ds_cfg: dict,
    vanilla_bin: Path,
    stream_bin: Path,
    threads: int,
    bench_root: Path,
    tmp_spill: Path,
    env: dict,
) -> dict:
    ds_id = ds_cfg["id"]
    print(f"\n=== Real Benchmark: {ds_id} ({ds_cfg['title']}) ===", flush=True)
    ddir = bench_root / ds_id
    if ddir.exists():
        shutil.rmtree(ddir)
    ddir.mkdir(parents=True)

    read_path = Path(ds_cfg["reads"])
    hic_paths = tuple(Path(x) for x in ds_cfg["hic"]) if ds_cfg.get("hic") else None
    staged_input_bytes = read_path.stat().st_size
    if hic_paths:
        staged_input_bytes += hic_paths[0].stat().st_size + hic_paths[1].stat().st_size

    # 1. Vanilla baseline (regular files on disk, default .bin writes)
    print(f"  [1/5] Running hifiasm_vanilla baseline (-t {threads}) ...", flush=True)
    vdir = ddir / "vanilla_disk"
    vprefix = vdir / "asm"
    v_args = ["-o", str(vprefix), "-t", str(threads)] + list(ds_cfg["flags"])
    if hic_paths:
        v_args += ["--h1", str(hic_paths[0]), "--h2", str(hic_paths[1])]
    v_args.append(str(read_path))
    v_wall, v_log, v_time = run_timed_vanilla(vanilla_bin, v_args, vdir, env)
    v_outs = collect_outputs(vprefix)
    v_disk = measure_disk_artifacts(vprefix)
    v_met = parse_log_metrics(v_log, v_time)
    shutil.rmtree(vdir)
    print(
        f"        vanilla: wall={v_wall:.2f}s pass1+2={v_met['pass12_sum_sec']}s "
        f"rss={v_met['peak_rss_gb']:.3f}GB .bin={v_disk['checkpoint_bin_gb']:.3f}GB ({v_disk['checkpoint_bin_count']} files)",
        flush=True,
    )

    def compare_outputs(cand_outs: dict, ref_outs: dict) -> dict:
        per_file = {}
        all_match = ( set(cand_outs.keys()) == set(ref_outs.keys()) )
        for k in sorted(set(ref_outs.keys()) | set(cand_outs.keys())):
            c_sha = cand_outs.get(k, {}).get("sha256")
            r_sha = ref_outs.get(k, {}).get("sha256")
            matched = (c_sha is not None and c_sha == r_sha)
            if not matched:
                all_match = False
            per_file[k] = {
                "matched": matched,
                "stream_sha256": c_sha,
                "vanilla_sha256": r_sha,
                "bytes": cand_outs.get(k, {}).get("bytes", 0),
            }
        return {"all_identical": all_match, "per_file": per_file}

    variants = {}

    # 2. Stream Pure-RAM (compressed input stream via stdin/FIFOs + --no-write-bin)
    print(f"  [2/5] Running stream_pure_ram_gz (compressed FIFO -> stdin, --stream-ram-gb 64 --no-write-bin) ...", flush=True)
    s1_dir = ddir / "stream_pure_ram_gz"
    s1_prefix = s1_dir / "asm"
    s1_args = ["-o", str(s1_prefix), "-t", str(threads), "--no-write-bin", "--stream-ram-gb", "64"] + list(ds_cfg["flags"])
    s1_wall, s1_log, s1_time = run_timed_stream(stream_bin, s1_args, read_path, s1_dir, env, decompress_upstream=False, hic_srcs=hic_paths)
    s1_outs = collect_outputs(s1_prefix)
    s1_disk = measure_disk_artifacts(s1_prefix)
    s1_met = parse_log_metrics(s1_log, s1_time)
    s1_cmp = compare_outputs(s1_outs, v_outs)
    assert s1_disk["checkpoint_bin_bytes"] == 0
    shutil.rmtree(s1_dir)
    variants["stream_pure_ram_gz"] = {
        "wall_sec": round(s1_wall, 3),
        "speedup_vs_vanilla": round(v_wall / s1_wall, 4),
        "metrics": s1_met,
        "disk": s1_disk,
        "sha256_exact_match": s1_cmp["all_identical"],
        "sha256_comparison": s1_cmp["per_file"],
    }
    print(
        f"        pure_ram_gz: wall={s1_wall:.2f}s pass1+2={s1_met['pass12_sum_sec']}s "
        f"rss={s1_met['peak_rss_gb']:.3f}GB store_ram={s1_met['stream_store']['peak_ram_gb']:.4f}GB SHA256_ALL={s1_cmp['all_identical']}",
        flush=True,
    )

    # 3. Stream Pure-RAM with live upstream decompression pipe (pigz -dc | hifiasm - --no-write-bin)
    print(f"  [3/5] Running stream_upstream_pipe (pigz -dc -> stdin/FIFOs, --stream-ram-gb 64 --no-write-bin) ...", flush=True)
    s2_dir = ddir / "stream_upstream_pipe"
    s2_prefix = s2_dir / "asm"
    s2_args = ["-o", str(s2_prefix), "-t", str(threads), "--no-write-bin", "--stream-ram-gb", "64"] + list(ds_cfg["flags"])
    s2_wall, s2_log, s2_time = run_timed_stream(stream_bin, s2_args, read_path, s2_dir, env, decompress_upstream=True, hic_srcs=hic_paths)
    s2_outs = collect_outputs(s2_prefix)
    s2_disk = measure_disk_artifacts(s2_prefix)
    s2_met = parse_log_metrics(s2_log, s2_time)
    s2_cmp = compare_outputs(s2_outs, v_outs)
    assert s2_disk["checkpoint_bin_bytes"] == 0
    shutil.rmtree(s2_dir)
    variants["stream_upstream_pipe"] = {
        "wall_sec": round(s2_wall, 3),
        "speedup_vs_vanilla": round(v_wall / s2_wall, 4),
        "metrics": s2_met,
        "disk": s2_disk,
        "sha256_exact_match": s2_cmp["all_identical"],
        "sha256_comparison": s2_cmp["per_file"],
    }
    print(
        f"        upstream_pipe: wall={s2_wall:.2f}s pass1+2={s2_met['pass12_sum_sec']}s "
        f"rss={s2_met['peak_rss_gb']:.3f}GB SHA256_ALL={s2_cmp['all_identical']}",
        flush=True,
    )

    # 4. Stream Pure O_TMPFILE Spool (--stream-ram-gb 0 --no-write-bin -> zero persistent StreamStore RAM!)
    print(f"  [4/5] Running stream_pure_spool (pigz -dc -> stdin/FIFOs, --stream-ram-gb 0 --no-write-bin) ...", flush=True)
    s3_dir = ddir / "stream_pure_spool"
    s3_prefix = s3_dir / "asm"
    s3_args = ["-o", str(s3_prefix), "-t", str(threads), "--no-write-bin", "--stream-ram-gb", "0"] + list(ds_cfg["flags"])
    s3_wall, s3_log, s3_time = run_timed_stream(stream_bin, s3_args, read_path, s3_dir, env, decompress_upstream=True, hic_srcs=hic_paths)
    s3_outs = collect_outputs(s3_prefix)
    s3_disk = measure_disk_artifacts(s3_prefix)
    s3_met = parse_log_metrics(s3_log, s3_time)
    s3_cmp = compare_outputs(s3_outs, v_outs)
    assert s3_met["stream_store"]["peak_ram_gb"] == 0.0
    assert s3_met["stream_store"]["spilled_disk_gb"] > 0.0
    assert len(list(tmp_spill.iterdir())) == 0, f"Leaked visible files in TMPDIR: {list(tmp_spill.iterdir())}"
    shutil.rmtree(s3_dir)
    variants["stream_pure_spool"] = {
        "wall_sec": round(s3_wall, 3),
        "speedup_vs_vanilla": round(v_wall / s3_wall, 4),
        "metrics": s3_met,
        "disk": s3_disk,
        "sha256_exact_match": s3_cmp["all_identical"],
        "sha256_comparison": s3_cmp["per_file"],
    }
    print(
        f"        pure_spool: wall={s3_wall:.2f}s pass1+2={s3_met['pass12_sum_sec']}s "
        f"rss={s3_met['peak_rss_gb']:.3f}GB store_ram=0.0000GB spilled={s3_met['stream_store']['spilled_disk_gb']:.4f}GB SHA256_ALL={s3_cmp['all_identical']}",
        flush=True,
    )

    # 5. Stream Hybrid RAM + O_TMPFILE Spool (~55% of StreamStore in RAM, ~45% spilled to O_TMPFILE)
    store_total_gb = s1_met["stream_store"]["peak_ram_gb"]
    hybrid_cap = f"{max(0.002, store_total_gb * 0.55):.5f}"
    print(f"  [5/5] Running stream_hybrid_spool (--stream-ram-gb {hybrid_cap} --no-write-bin) ...", flush=True)
    s4_dir = ddir / "stream_hybrid_spool"
    s4_prefix = s4_dir / "asm"
    s4_args = ["-o", str(s4_prefix), "-t", str(threads), "--no-write-bin", "--stream-ram-gb", hybrid_cap] + list(ds_cfg["flags"])
    s4_wall, s4_log, s4_time = run_timed_stream(stream_bin, s4_args, read_path, s4_dir, env, decompress_upstream=True, hic_srcs=hic_paths)
    s4_outs = collect_outputs(s4_prefix)
    s4_disk = measure_disk_artifacts(s4_prefix)
    s4_met = parse_log_metrics(s4_log, s4_time)
    s4_cmp = compare_outputs(s4_outs, v_outs)
    if store_total_gb > 0.015:
        assert s4_met["stream_store"]["peak_ram_gb"] > 0.0 and s4_met["stream_store"]["spilled_disk_gb"] > 0.0
    assert len(list(tmp_spill.iterdir())) == 0
    shutil.rmtree(s4_dir)
    variants["stream_hybrid_spool"] = {
        "hybrid_cap_gb": float(hybrid_cap),
        "wall_sec": round(s4_wall, 3),
        "speedup_vs_vanilla": round(v_wall / s4_wall, 4),
        "metrics": s4_met,
        "disk": s4_disk,
        "sha256_exact_match": s4_cmp["all_identical"],
        "sha256_comparison": s4_cmp["per_file"],
    }
    print(
        f"        hybrid_spool: wall={s4_wall:.2f}s rss={s4_met['peak_rss_gb']:.3f}GB "
        f"store_ram={s4_met['stream_store']['peak_ram_gb']:.4f}GB spilled={s4_met['stream_store']['spilled_disk_gb']:.4f}GB SHA256_ALL={s4_cmp['all_identical']}",
        flush=True,
    )

    shutil.rmtree(ddir)
    return {
        "id": ds_id,
        "title": ds_cfg["title"],
        "threads": threads,
        "staged_input_bytes": staged_input_bytes,
        "staged_input_gb": round(staged_input_bytes / (1024 ** 3), 4),
        "vanilla": {
            "wall_sec": round(v_wall, 3),
            "metrics": v_met,
            "disk": v_disk,
        },
        "variants": variants,
        "verified_outputs": v_outs,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vanilla", required=True, type=Path)
    ap.add_argument("--stream-bin", required=True, type=Path)
    ap.add_argument("--data-dir", required=True, type=Path)
    ap.add_argument("--tmpdir", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=24)
    ap.add_argument("--out-json", required=True, type=Path)
    args = ap.parse_args()

    args.tmpdir.mkdir(parents=True, exist_ok=True)
    bench_root = args.tmpdir / f"real_bench_{os.getpid()}"
    if bench_root.exists():
        shutil.rmtree(bench_root)
    bench_root.mkdir(parents=True)
    tmp_spill = bench_root / "spill_tmp"
    tmp_spill.mkdir(parents=True)

    env = os.environ.copy()
    env["TMPDIR"] = str(tmp_spill)

    datasets = [
        {
            "id": "human_chr11_2m_hifi",
            "title": "Human chr11 2 Mbp PacBio HiFi (chr11-2M.fa.gz)",
            "reads": args.data_dir / "chr11-2M.fa.gz",
            "hic": None,
            "flags": ["-f0"],
        },
        {
            "id": "dmel_ont_r1041_25x",
            "title": "Drosophila melanogaster idDroMela1 Whole-Genome ~23x ONT R10.4.1 SUP Simplex (--ont)",
            "reads": args.data_dir / "dmel_ont_r1041_25x.fq.gz",
            "hic": None,
            "flags": ["--ont"],
        },
        {
            "id": "dmel_hifi_30x_hic",
            "title": "Drosophila melanogaster idDroMela1 Whole-Genome ~30x PacBio HiFi + 2M Dovetail Hi-C Pairs",
            "reads": args.data_dir / "dmel_hifi_30x.fq.gz",
            "hic": (args.data_dir / "dmel_hic_R1.fq.gz", args.data_dir / "dmel_hic_R2.fq.gz"),
            "flags": [],
        },
    ]

    all_results = []
    for ds_cfg in datasets:
        res = run_dataset_benchmark(
            ds_cfg,
            args.vanilla,
            args.stream_bin,
            args.threads,
            bench_root,
            tmp_spill,
            env,
        )
        all_results.append(res)
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps({"benchmarks": all_results}, indent=2))

    shutil.rmtree(bench_root)
    print(f"\n[SUCCESS] Completed all {len(all_results)} real-genome benchmark suites with 100% SHA-256 bit-identity!", flush=True)


if __name__ == "__main__":
    main()
