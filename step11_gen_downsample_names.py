#!/usr/bin/env python3
"""step11：双端配对 read-name 列表；shuffle 后嵌套切片（低深度 ⊂ 高深度）。"""

from __future__ import annotations

import argparse
import hashlib
import os
import random
import subprocess
import sys
import tempfile
from typing import List, Tuple


def depth_label(d: float) -> str:
    return f"{d:.1f}"


def load_paired_ids(r1: str, r2: str, tmp_path: str) -> Tuple[List[str], int]:
    """shell paste+awk 流式提取配对 base ID（比 Python 逐条读 fq 更快）。"""
    script = f"""
set -euo pipefail
AWK='NR%4==1{{line=$0; sub(/^@/, "", line); split(line, a, " "); id=a[1]; sub(/\\/[12]$/, "", id); print id}}'
paste \\
  <( (pigz -dc "{r1}" 2>/dev/null || gzip -dc "{r1}") | awk "$AWK" ) \\
  <( (pigz -dc "{r2}" 2>/dev/null || gzip -dc "{r2}") | awk "$AWK" ) \\
| awk -F'\\t' 'BEGIN{{mis=0}} {{
    if ($1 != $2) {{ mis++; next }}
    print $1
}} END {{
    if (mis > 0) {{ print "pair_mismatch=" mis > "/dev/stderr"; exit 3 }}
}}' > "{tmp_path}"
"""
    proc = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
    if proc.returncode == 3:
        mis = proc.stderr.strip()
        raise RuntimeError(f"R1/R2 name mismatch: {mis}")
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"load paired ids failed rc={proc.returncode}")

    ids: List[str] = []
    with open(tmp_path, "rt") as fh:
        for line in fh:
            s = line.strip()
            if s:
                ids.append(s)
    return ids, 0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--names-dir", required=True)
    ap.add_argument("--expected-depth", type=float, required=True)
    ap.add_argument("--depths", required=True)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    if not os.path.isfile(args.r1):
        print(f"ERROR: missing R1 {args.r1}", file=sys.stderr)
        sys.exit(2)
    if not os.path.isfile(args.r2):
        print(f"ERROR: missing R2 {args.r2}", file=sys.stderr)
        sys.exit(2)

    depths = sorted({float(x) for x in args.depths.split()})
    out_dir = os.path.join(args.names_dir, args.sample)
    os.makedirs(out_dir, exist_ok=True)

    tmp_ids = os.path.join(out_dir, f".{args.sample}.paired_ids.tmp")
    try:
        print(f"[{args.sample}] loading paired ids ...", file=sys.stderr)
        ids, _ = load_paired_ids(args.r1, args.r2, tmp_ids)
    finally:
        if os.path.isfile(tmp_ids):
            os.remove(tmp_ids)

    n_pairs = len(ids)
    if n_pairs <= 0:
        print(f"ERROR: no paired reads in {args.r1}", file=sys.stderr)
        sys.exit(2)

    seed = args.seed + int(hashlib.md5(args.sample.encode()).hexdigest()[:8], 16)
    rng = random.Random(seed)
    rng.shuffle(ids)

    ks = {}
    for depth in depths:
        frac = depth / args.expected_depth
        ks[depth] = max(1, min(n_pairs, round(n_pairs * frac)))

    prev_k = 0
    for depth in depths:
        k = ks[depth]
        if k < prev_k:
            print(f"ERROR: non-nested k depth={depth} k={k} prev={prev_k}", file=sys.stderr)
            sys.exit(2)
        prev_k = k

    for depth in depths:
        k = ks[depth]
        dstr = depth_label(depth)
        out_path = os.path.join(out_dir, f"names_{dstr}x.txt")
        subset = ids[:k]
        with open(out_path, "w") as fh:
            fh.write("\n".join(subset))
            fh.write("\n")
        print(
            f"[{args.sample}] names_{dstr}x n={len(subset)}/{n_pairs} "
            f"target={depth}x nested=1 seed={seed} -> {out_path}",
            file=sys.stderr,
        )

    sys.exit(0)


if __name__ == "__main__":
    main()
