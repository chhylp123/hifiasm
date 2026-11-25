#!/usr/bin/env python3
"""
Automated ASan/UBSan + SHA-256 bit-identity verification suite for single-pass streaming hifiasm.
Verifies ALL 5 assembly modes (HiFi -l0, HiFi --primary, Trio -1/-2, Hi-C --h1/--h2 FIFOs, ONT --ont)
at BOTH -t 1 and -t 16 across ALL 3 StreamStore policies:
  1. Pure RAM mode (--stream-ram-gb 64)
  2. Pure O_TMPFILE spool mode (--stream-ram-gb 0, zero persistent StreamStore RAM)
  3. Hybrid RAM + O_TMPFILE spool mode (--stream-ram-gb <batch1_cap>, multi-batch split across RAM and TMPDIR)
plus -fsanitize=address,undefined (ASan + LSan + UBSan) against frozen hifiasm_vanilla.
"""

import argparse
import gzip
import hashlib
import json
import os
import random
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


def revcomp(seq: str) -> str:
    tr = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tr)[::-1]


def generate_diploid_dataset(data_dir: Path, make_yak_bin: Path) -> dict:
    data_dir.mkdir(parents=True, exist_ok=True)
    stamp = data_dir / "dataset_ready_v5.json"
    if stamp.exists():
        return json.loads(stamp.read_text())

    rng = random.Random(20261124)
    bases = "ACGT"
    # Single 720 kb diploid chromosome (1.44 Mb diploid -> ~25.9 Mbp HiFi & ~31.7 Mbp ONT = 2 batches)
    # with two asymmetric homozygous intervals so superbubbles form cleanly with distinct link weights.
    chr_lens = [720000]
    hap1_chrs = []
    hap2_chrs = []

    for c_idx, clen in enumerate(chr_lens):
        h1_list = [bases[rng.randrange(4)] for _ in range(clen)]
        if c_idx == 0:
            for pos in range(45000, 45012):
                h1_list[pos] = "N"
        h1 = "".join(h1_list)

        h2_parts = []
        i = 0
        while i < len(h1):
            ch = h1[i]
            in_hom = (180000 <= i < 199000) or (490000 <= i < 514000)
            if ch == "N" or in_hom:
                h2_parts.append(ch)
                i += 1
                continue
            r = rng.random()
            if r < 0.010:
                alts = [b for b in bases if b != ch]
                h2_parts.append(alts[rng.randrange(3)])
                i += 1
            elif r < 0.011 and i + 5 < len(h1):
                for _ in range(3):
                    h2_parts.append(bases[rng.randrange(4)])
                h2_parts.append(ch)
                i += 1
            elif r < 0.012 and i + 5 < len(h1):
                i += 2
            else:
                h2_parts.append(ch)
                i += 1
        hap1_chrs.append(h1)
        hap2_chrs.append("".join(h2_parts))

    # 1. PacBio HiFi FASTQ.gz (~36x diploid = ~25.9 Mbp -> 2 kt_pipeline batches!)
    hifi_fq = data_dir / "synth_hifi.fq.gz"
    n_hifi = 0
    total_hifi_bp = 0
    with gzip.open(hifi_fq, "wt", compresslevel=1) as fh:
        for hap_idx, chrs in enumerate([hap1_chrs, hap2_chrs]):
            for c_idx, hseq in enumerate(chrs):
                target_bp = len(hseq) * 18
                acc_bp = 0
                while acc_bp < target_bp:
                    rlen = rng.randint(11500, 15500)
                    start = rng.randint(0, len(hseq) - rlen)
                    sub = hseq[start : start + rlen]
                    if rng.random() < 0.5:
                        sub = revcomp(sub)
                    s_arr = list(sub)
                    q_arr = []
                    for idx_b, b in enumerate(s_arr):
                        if b == "N":
                            q_arr.append("!")
                        elif rng.random() < 0.0015:
                            s_arr[idx_b] = bases[rng.randrange(4)]
                            q_arr.append(chr(33 + rng.randint(10, 16)))
                        else:
                            q_arr.append(chr(33 + rng.randint(28, 40)))
                    if n_hifi % 25 == 0 and len(s_arr) > 20:
                        s_arr[10] = s_arr[10].lower()
                    fh.write(f"@hifi_h{hap_idx}_c{c_idx}_{n_hifi} length={rlen}\n")
                    fh.write("".join(s_arr) + "\n+\n" + "".join(q_arr) + "\n")
                    n_hifi += 1
                    acc_bp += rlen
                    total_hifi_bp += rlen

    # 2. ONT R10.4 Simplex FASTQ.gz (~44x diploid = ~31.7 Mbp -> 2 kt_pipeline batches!)
    ont_fq = data_dir / "synth_ont.fq.gz"
    n_ont = 0
    total_ont_bp = 0
    with gzip.open(ont_fq, "wt", compresslevel=1) as fh:
        for k in range(8):
            bad_s = "".join(bases[rng.randrange(4)] for _ in range(600 if k < 4 else 2500))
            bad_q = chr(33 + 5) * len(bad_s)
            fh.write(f"@ont_filtered_{k}\n{bad_s}\n+\n{bad_q}\n")

        for hap_idx, chrs in enumerate([hap1_chrs, hap2_chrs]):
            for c_idx, hseq in enumerate(chrs):
                target_bp = len(hseq) * 22
                acc_bp = 0
                while acc_bp < target_bp:
                    rlen = rng.randint(10000, 16000)
                    start = rng.randint(0, len(hseq) - rlen)
                    sub = hseq[start : start + rlen]
                    if rng.random() < 0.5:
                        sub = revcomp(sub)
                    s_arr = []
                    q_arr = []
                    for b in sub:
                        if b == "N":
                            s_arr.append("N")
                            q_arr.append("!")
                            continue
                        er = rng.random()
                        if er < 0.014:
                            alts = [x for x in bases if x != b]
                            s_arr.append(alts[rng.randrange(3)])
                            q_arr.append(chr(33 + rng.randint(5, 9)))
                        elif er < 0.018:
                            s_arr.append(bases[rng.randrange(4)])
                            q_arr.append(chr(33 + rng.randint(5, 8)))
                            s_arr.append(b)
                            q_arr.append(chr(33 + rng.randint(18, 30)))
                        elif er < 0.022:
                            continue
                        else:
                            s_arr.append(b)
                            q_arr.append(chr(33 + rng.randint(16, 36)))
                    fh.write(f"@ont_h{hap_idx}_c{c_idx}_{n_ont}\n")
                    fh.write("".join(s_arr) + "\n+\n" + "".join(q_arr) + "\n")
                    n_ont += 1
                    acc_bp += len(s_arr)
                    total_ont_bp += len(s_arr)

    # 3. Parental k-mer tables (pat.yak, mat.yak)
    pat_fa = data_dir / "pat_reads.fa.gz"
    mat_fa = data_dir / "mat_reads.fa.gz"
    for out_fa, chrs, tag in [(pat_fa, hap1_chrs, "pat"), (mat_fa, hap2_chrs, "mat")]:
        with gzip.open(out_fa, "wt", compresslevel=1) as fh:
            rid = 0
            for hseq in chrs:
                for _ in range(24):
                    for pos in range(0, len(hseq) - 150, 75):
                        sub = hseq[pos : pos + 150]
                        fh.write(f">{tag}_{rid}\n{sub}\n")
                        rid += 1
    pat_yak = data_dir / "pat.yak"
    mat_yak = data_dir / "mat.yak"
    subprocess.run([str(make_yak_bin), str(pat_yak), str(pat_fa)], check=True)
    subprocess.run([str(make_yak_bin), str(mat_yak), str(mat_fa)], check=True)

    # 4. Paired-end Hi-C reads (88% cis intra-haplotype + 12% trans-haplotype noise)
    all_hap_seqs = []
    for hap_idx, chrs in enumerate([hap1_chrs, hap2_chrs]):
        for c_idx, hseq in enumerate(chrs):
            all_hap_seqs.append((hap_idx, c_idx, hseq))
    hic1_fq = data_dir / "synth_hic_R1.fq.gz"
    hic2_fq = data_dir / "synth_hic_R2.fq.gz"
    n_hic = 0
    pairs_per_seq = 32000 // len(all_hap_seqs)
    with gzip.open(hic1_fq, "wt", compresslevel=1) as f1, gzip.open(hic2_fq, "wt", compresslevel=1) as f2:
        for hap_idx, c_idx, hseq in all_hap_seqs:
            for _ in range(pairs_per_seq):
                if rng.random() < 0.88:
                    span = rng.randint(1500, 320000)
                    p1 = rng.randint(0, len(hseq) - span - 160)
                    p2 = p1 + span
                    r1 = hseq[p1 : p1 + 150]
                    r2 = revcomp(hseq[p2 : p2 + 150])
                else:
                    other_hap = 1 - hap_idx
                    other_seq = hap2_chrs[c_idx] if other_hap == 1 else hap1_chrs[c_idx]
                    p1 = rng.randint(0, len(hseq) - 160)
                    p2 = rng.randint(0, len(other_seq) - 160)
                    r1 = hseq[p1 : p1 + 150]
                    r2 = revcomp(other_seq[p2 : p2 + 150])
                if rng.random() < 0.5:
                    r1, r2 = r2, r1
                q = "I" * 150
                f1.write(f"@hic_h{hap_idx}_c{c_idx}_{n_hic}/1\n{r1}\n+\n{q}\n")
                f2.write(f"@hic_h{hap_idx}_c{c_idx}_{n_hic}/2\n{r2}\n+\n{q}\n")
                n_hic += 1

    meta = {
        "hap1_len": sum(len(x) for x in hap1_chrs),
        "hap2_len": sum(len(x) for x in hap2_chrs),
        "hifi_fq": str(hifi_fq),
        "ont_fq": str(ont_fq),
        "pat_yak": str(pat_yak),
        "mat_yak": str(mat_yak),
        "hic1_fq": str(hic1_fq),
        "hic2_fq": str(hic2_fq),
        "n_hifi_reads": n_hifi,
        "total_hifi_bp": total_hifi_bp,
        "n_ont_reads": n_ont,
        "total_ont_bp": total_ont_bp,
        "n_hic_pairs": n_hic,
    }
    stamp.write_text(json.dumps(meta, indent=2))
    return meta


def collect_output_hashes(prefix: Path) -> dict:
    parent = prefix.parent
    base = prefix.name
    hashes = {}
    for p in sorted(parent.glob(f"{base}.*")):
        if p.name.endswith(".gfa") or p.name.endswith(".bed"):
            hashes[p.name[len(base) + 1 :]] = {
                "sha256": sha256_file(p),
                "bytes": p.stat().st_size,
            }
    return hashes


def count_bin_files(prefix: Path) -> list:
    return [p.name for p in sorted(prefix.parent.glob(f"{prefix.name}.*.bin"))]


def parse_stream_metrics(log_text: str) -> dict:
    m = re.search(
        r"StreamStore captured (\d+) reads \((\d+) bp\) peak_ram=([0-9.]+) GB spilled_disk=([0-9.]+) GB",
        log_text,
    )
    if not m:
        raise RuntimeError(f"Could not find StreamStore log line in output:\n{log_text[-1500:]}")
    return {
        "reads": int(m.group(1)),
        "bp": int(m.group(2)),
        "peak_ram_gb": float(m.group(3)),
        "spilled_disk_gb": float(m.group(4)),
    }


def run_streaming_fifo_cmd(
    exe: Path,
    args_before_reads: list,
    read_src: Path,
    work_dir: Path,
    env: dict,
    use_stdin: bool = False,
    hic_srcs: tuple = None,
) -> tuple:
    work_dir.mkdir(parents=True, exist_ok=True)
    writers = []
    cmd = [str(exe)] + list(args_before_reads)

    if hic_srcs is not None:
        fifo_h1 = work_dir / "in_h1.fifo"
        fifo_h2 = work_dir / "in_h2.fifo"
        os.mkfifo(fifo_h1)
        os.mkfifo(fifo_h2)
        w1 = subprocess.Popen(["bash", "-c", f'cat "{hic_srcs[0]}" > "{fifo_h1}"'], env=env)
        w2 = subprocess.Popen(["bash", "-c", f'cat "{hic_srcs[1]}" > "{fifo_h2}"'], env=env)
        writers.extend([w1, w2])
        cmd.extend(["--h1", str(fifo_h1), "--h2", str(fifo_h2)])

    t0 = time.perf_counter()
    if use_stdin:
        cmd.append("-")
        log_path = work_dir / "run.log"
        with open(read_src, "rb") as fin, open(log_path, "wb") as flog:
            proc = subprocess.Popen(cmd, stdin=fin, stdout=flog, stderr=subprocess.STDOUT, env=env)
            rc = proc.wait()
    else:
        fifo_r = work_dir / "in_reads.fifo"
        os.mkfifo(fifo_r)
        wr = subprocess.Popen(["bash", "-c", f'cat "{read_src}" > "{fifo_r}"'], env=env)
        writers.append(wr)
        cmd.append(str(fifo_r))
        log_path = work_dir / "run.log"
        with open(log_path, "wb") as flog:
            proc = subprocess.Popen(cmd, stdout=flog, stderr=subprocess.STDOUT, env=env)
            rc = proc.wait()

    for w in writers:
        wrc = w.wait()
        if wrc != 0:
            raise RuntimeError(f"FIFO writer exited with non-zero status {wrc} in {work_dir}")
    elapsed = time.perf_counter() - t0
    if rc != 0:
        tail = (work_dir / "run.log").read_text(errors="replace")[-3000:]
        raise RuntimeError(f"Command failed (rc={rc}): {' '.join(cmd)}\n{tail}")
    return elapsed, (work_dir / "run.log").read_text(errors="replace")


def run_vanilla_cmd(exe: Path, args: list, work_dir: Path, env: dict) -> tuple:
    work_dir.mkdir(parents=True, exist_ok=True)
    log_path = work_dir / "run.log"
    t0 = time.perf_counter()
    with open(log_path, "wb") as flog:
        proc = subprocess.run([str(exe)] + args, stdout=flog, stderr=subprocess.STDOUT, env=env)
    elapsed = time.perf_counter() - t0
    if proc.returncode != 0:
        tail = log_path.read_text(errors="replace")[-3000:]
        raise RuntimeError(f"Vanilla failed (rc={proc.returncode}):\n{tail}")
    return elapsed, log_path.read_text(errors="replace")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vanilla", required=True, type=Path)
    ap.add_argument("--stream-bin", required=True, type=Path)
    ap.add_argument("--san-bin", required=True, type=Path)
    ap.add_argument("--make-yak", required=True, type=Path)
    ap.add_argument("--tmpdir", required=True, type=Path)
    ap.add_argument("--out-json", required=True, type=Path)
    args = ap.parse_args()

    args.tmpdir.mkdir(parents=True, exist_ok=True)
    suite_root = args.tmpdir / f"suite_run_{os.getpid()}"
    if suite_root.exists():
        shutil.rmtree(suite_root)
    suite_root.mkdir(parents=True)
    tmp_spill = suite_root / "spill_tmp"
    tmp_spill.mkdir(parents=True)

    base_env = os.environ.copy()
    base_env["TMPDIR"] = str(tmp_spill)

    san_env = base_env.copy()
    san_env["ASAN_OPTIONS"] = "detect_leaks=1:halt_on_error=1:abort_on_error=1"
    san_env["UBSAN_OPTIONS"] = "print_stacktrace=1:halt_on_error=1:abort_on_error=1"

    data_dir = args.tmpdir / "synth_diploid_v2"
    ds = generate_diploid_dataset(data_dir, args.make_yak)
    print(f"[INFO] Generated multi-batch diploid dataset: {json.dumps(ds)}", flush=True)

    # Batch 0 (~20 Mbp) compresses to ~10-12 MB for HiFi/ONT FASTQ; cap=0.015 GiB (~16.1 MB) keeps Batch 0 in RAM and spills Batch 1 to O_TMPFILE!
    modes = [
        {
            "id": "hifi_l0",
            "label": "PacBio HiFi haploid (-l0)",
            "common_flags": ["-f22", "-l0"],
            "read_file": Path(ds["hifi_fq"]),
            "hic": None,
            "hybrid_cap_gb": "0.015",
        },
        {
            "id": "hifi_primary",
            "label": "PacBio HiFi diploid (--primary)",
            "common_flags": ["-f22", "--primary"],
            "read_file": Path(ds["hifi_fq"]),
            "hic": None,
            "hybrid_cap_gb": "0.015",
        },
        {
            "id": "trio_yak",
            "label": "Trio diploid (-1 pat.yak -2 mat.yak)",
            "common_flags": ["-f22", "-1", ds["pat_yak"], "-2", ds["mat_yak"]],
            "read_file": Path(ds["hifi_fq"]),
            "hic": None,
            "hybrid_cap_gb": "0.015",
        },
        {
            "id": "hic_fifo",
            "label": "Hi-C phased (--h1/--h2 FIFOs)",
            "common_flags": ["-f22"],
            "read_file": Path(ds["hifi_fq"]),
            "hic": (Path(ds["hic1_fq"]), Path(ds["hic2_fq"])),
            "hybrid_cap_gb": "0.015",
        },
        {
            "id": "ont_simplex",
            "label": "ONT R10.4 Simplex (--ont with base qualities)",
            "common_flags": ["-f22", "--ont", "--hg-size", "700k"],
            "read_file": Path(ds["ont_fq"]),
            "hic": None,
            "hybrid_cap_gb": "0.018",
        },
    ]

    results = []
    for mode in modes:
        for threads in [1, 16]:
            case_key = f"{mode['id']}_t{threads}"
            print(f"[RUN] {case_key}: {mode['label']} (-t {threads}) ...", flush=True)
            cdir = suite_root / case_key

            # 1. Frozen hifiasm_vanilla on regular files
            vdir = cdir / "vanilla"
            vprefix = vdir / "asm"
            v_args = ["-o", str(vprefix), "-t", str(threads)] + list(mode["common_flags"])
            if mode["hic"]:
                v_args += ["--h1", str(mode["hic"][0]), "--h2", str(mode["hic"][1])]
            v_args.append(str(mode["read_file"]))
            v_sec, _ = run_vanilla_cmd(args.vanilla, v_args, vdir, base_env)
            v_hashes = collect_output_hashes(vprefix)
            v_bins = count_bin_files(vprefix)
            assert len(v_hashes) > 0, f"No GFA/BED outputs produced for {case_key}"
            assert any(
                v["bytes"] > 1000
                for k, v in v_hashes.items()
                if k == "p_ctg.gfa" or k.endswith(".p_ctg.gfa")
            ), f"Empty primary contigs for {case_key}: {v_hashes}"
            assert len(v_bins) > 0, f"Expected vanilla .bin files for {case_key}"

            # 2. Policy A: Pure RAM StreamStore (--stream-ram-gb 64) via stdin (-) or FIFO + --no-write-bin
            sdir = cdir / "stream_ram"
            sprefix = sdir / "asm"
            s_args = [
                "-o", str(sprefix),
                "-t", str(threads),
                "--no-write-bin",
                "--stream-ram-gb", "64",
            ] + list(mode["common_flags"])
            use_stdin = (threads == 1 and mode["hic"] is None)
            s_sec, s_log = run_streaming_fifo_cmd(
                args.stream_bin,
                s_args,
                mode["read_file"],
                sdir,
                base_env,
                use_stdin=use_stdin,
                hic_srcs=mode["hic"],
            )
            s_hashes = collect_output_hashes(sprefix)
            s_bins = count_bin_files(sprefix)
            s_met = parse_stream_metrics(s_log)
            assert s_hashes == v_hashes, f"SHA-256 mismatch in Pure RAM mode for {case_key}"
            assert len(s_bins) == 0, f"Expected 0 .bin files with --no-write-bin, got {s_bins}"
            assert s_met["peak_ram_gb"] > 0.0 and s_met["spilled_disk_gb"] == 0.0, f"Unexpected RAM metrics: {s_met}"

            # 3. Policy B: Pure O_TMPFILE Spool (--stream-ram-gb 0 -> 0 B persistent StreamStore RAM!)
            pdir = cdir / "stream_pure_spool"
            pprefix = pdir / "asm"
            p_args = [
                "-o", str(pprefix),
                "-t", str(threads),
                "--no-write-bin",
                "--stream-ram-gb", "0",
            ] + list(mode["common_flags"])
            p_sec, p_log = run_streaming_fifo_cmd(
                args.stream_bin,
                p_args,
                mode["read_file"],
                pdir,
                base_env,
                use_stdin=(not use_stdin and mode["hic"] is None),
                hic_srcs=mode["hic"],
            )
            p_hashes = collect_output_hashes(pprefix)
            p_bins = count_bin_files(pprefix)
            p_met = parse_stream_metrics(p_log)
            assert p_hashes == v_hashes, f"SHA-256 mismatch in Pure O_TMPFILE spool mode (--stream-ram-gb 0) for {case_key}"
            assert len(p_bins) == 0
            assert p_met["peak_ram_gb"] == 0.0 and p_met["spilled_disk_gb"] > 0.0, (
                f"Expected peak_ram_gb == 0.0 and spilled_disk_gb > 0.0 for --stream-ram-gb 0, got {p_met}"
            )
            assert len(list(tmp_spill.iterdir())) == 0, f"Visible files leaked in TMPDIR: {list(tmp_spill.iterdir())}"

            # 4. Policy C: Hybrid RAM + O_TMPFILE Spool (Batch 0 in RAM, Batch 1+ in unlinked O_TMPFILE!)
            hybrid_cap_str = f"{s_met['peak_ram_gb'] * 0.85:.5f}"
            hdir = cdir / "stream_hybrid_spool"
            hprefix = hdir / "asm"
            h_args = [
                "-o", str(hprefix),
                "-t", str(threads),
                "--no-write-bin",
                "--stream-ram-gb", hybrid_cap_str,
            ] + list(mode["common_flags"])
            h_sec, h_log = run_streaming_fifo_cmd(
                args.stream_bin,
                h_args,
                mode["read_file"],
                hdir,
                base_env,
                use_stdin=False,
                hic_srcs=mode["hic"],
            )
            h_hashes = collect_output_hashes(hprefix)
            h_met = parse_stream_metrics(h_log)
            assert h_hashes == v_hashes, f"SHA-256 mismatch in Hybrid RAM+O_TMPFILE mode for {case_key}"
            assert h_met["peak_ram_gb"] > 0.0 and h_met["spilled_disk_gb"] > 0.0, (
                f"Expected BOTH peak_ram_gb > 0 and spilled_disk_gb > 0 in Hybrid mode, got {h_met}"
            )
            assert len(list(tmp_spill.iterdir())) == 0, f"Visible files leaked in TMPDIR: {list(tmp_spill.iterdir())}"

            # 5. Full ASan + LSan + UBSan verification across all 5 modes (exercises Hybrid RAM+O_TMPFILE + FIFOs)
            san_checked = False
            if threads == 16:
                adir = cdir / "san"
                aprefix = adir / "asm"
                a_args = [
                    "-o", str(aprefix),
                    "-t", "8",
                    "--no-write-bin",
                    "--stream-ram-gb", hybrid_cap_str,
                ] + list(mode["common_flags"])
                _, a_log = run_streaming_fifo_cmd(
                    args.san_bin,
                    a_args,
                    mode["read_file"],
                    adir,
                    san_env,
                    use_stdin=False,
                    hic_srcs=mode["hic"],
                )
                a_hashes = collect_output_hashes(aprefix)
                assert a_hashes == v_hashes, f"SHA-256 mismatch in ASan/UBSan run for {case_key}"
                assert "ERROR: AddressSanitizer" not in a_log
                assert "ERROR: LeakSanitizer" not in a_log
                assert "runtime error:" not in a_log
                san_checked = True

            shutil.rmtree(cdir)
            results.append({
                "case": case_key,
                "mode": mode["id"],
                "label": mode["label"],
                "threads": threads,
                "vanilla_sec": round(v_sec, 3),
                "stream_ram_sec": round(s_sec, 3),
                "stream_pure_spool_sec": round(p_sec, 3),
                "stream_hybrid_spool_sec": round(h_sec, 3),
                "metrics_pure_ram": s_met,
                "metrics_pure_spool": p_met,
                "metrics_hybrid_spool": h_met,
                "num_outputs_verified": len(v_hashes),
                "outputs": v_hashes,
                "vanilla_bin_count": len(v_bins),
                "stream_bin_count": len(s_bins),
                "asan_ubsan_verified": san_checked,
            })
            print(
                f"[PASS] {case_key}: {len(v_hashes)} GFA/BED files 100% SHA-256 identical across "
                f"Pure-RAM ({s_met['peak_ram_gb']:.3f}G/0G), "
                f"Pure-Spool (0G/{p_met['spilled_disk_gb']:.3f}G), "
                f"Hybrid ({h_met['peak_ram_gb']:.3f}G/{h_met['spilled_disk_gb']:.3f}G), "
                f"ASan/UBSan={san_checked}",
                flush=True,
            )

    shutil.rmtree(suite_root)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps({"dataset": ds, "results": results}, indent=2))
    print(f"[SUCCESS] All {len(results)} mode x thread configurations passed Pure-RAM + Pure-Spool + Hybrid-Spool SHA-256 bit-identity!", flush=True)


if __name__ == "__main__":
    main()
