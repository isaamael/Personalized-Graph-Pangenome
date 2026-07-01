#!/usr/bin/env python3
"""step03: 按 fragments 双坐标从亲本 FASTA 切片，拼接单倍型 FASTA。"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from functools import partial
from multiprocessing import Pool
from typing import Dict, List, Tuple

import pysam


class FastaCache:
    def __init__(self, mm_fa: str, ts_fa: str):
        self.mm = pysam.FastaFile(mm_fa)
        self.ts = pysam.FastaFile(ts_fa)

    def fetch(self, parent: str, chrom: str, start: int, end: int) -> str:
        if chrom in ("", "-") or start < 1 or end < start:
            return ""
        fa = self.mm if parent == "MM" else self.ts
        try:
            return fa.fetch(chrom, start - 1, end)
        except (KeyError, ValueError):
            return ""

    def close(self) -> None:
        self.mm.close()
        self.ts.close()


def extract_fragment(fa: FastaCache, frag: dict) -> str:
    if frag["parent"] == "MM":
        return fa.fetch("MM", frag["mm_chr"], frag["mm_lo"], frag["mm_hi"])
    return fa.fetch("TS", frag["ts_chr"], frag["ts_lo"], frag["ts_hi"])


def load_fragments(path: str) -> Tuple[str, str, List[dict]]:
    frags: List[dict] = []
    sample = gamete = ""
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            sample = row["sample"]
            gamete = row["gamete"]
            frags.append({
                "frag_idx": int(row["frag_idx"]),
                "parent": row["parent"],
                "mm_chr": row["mm_chr"],
                "mm_lo": int(row["mm_lo"]),
                "mm_hi": int(row["mm_hi"]),
                "ts_chr": row["ts_chr"],
                "ts_lo": int(row["ts_lo"]),
                "ts_hi": int(row["ts_hi"]),
            })
    frags.sort(key=lambda x: x["frag_idx"])
    return sample, gamete, frags


def build_haplotype_fasta(mm_fa: str, ts_fa: str, frags: List[dict], out_fa: str) -> None:
    fa = FastaCache(mm_fa, ts_fa)
    try:
        by_chr: Dict[str, List[dict]] = {}
        for f in frags:
            by_chr.setdefault(f["mm_chr"], []).append(f)
        for c in by_chr:
            by_chr[c].sort(key=lambda x: x["frag_idx"])

        with open(out_fa, "w") as out:
            for mm_chr in sorted(by_chr.keys(), key=lambda c: int(c.replace("chr", ""))):
                parts: List[str] = []
                for f in by_chr[mm_chr]:
                    seq = extract_fragment(fa, f)
                    if seq:
                        parts.append(seq)
                if parts:
                    out.write(f">{mm_chr}\n")
                    seq = "".join(parts)
                    for i in range(0, len(seq), 80):
                        out.write(seq[i:i + 80] + "\n")
    finally:
        fa.close()


def process_one_frag(fn: str, args: argparse.Namespace) -> str:
    path = os.path.join(args.truth_dir, fn)
    sample, gamete, frags = load_fragments(path)
    out_fa = os.path.join(args.out_dir, f"{sample}.gamete{gamete}.fa")
    if os.path.isfile(out_fa) and os.path.getsize(out_fa) > 0 and not args.force:
        return f"skip {out_fa}"
    build_haplotype_fasta(args.mm_fa, args.ts_fa, frags, out_fa)
    return f"wrote {out_fa}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-dir", required=True)
    ap.add_argument("--panel-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--mm-fa", required=True)
    ap.add_argument("--ts-fa", required=True)
    ap.add_argument("--sample", default="")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    frag_files = sorted(
        f for f in os.listdir(args.truth_dir)
        if f.endswith(".fragments.tsv")
        and (not args.sample or f.startswith(f"{args.sample}."))
    )
    if not frag_files:
        sys.exit("no fragments.tsv found")

    workers = max(1, min(args.workers, len(frag_files)))
    if workers == 1:
        for fn in frag_files:
            print(process_one_frag(fn, args), file=sys.stderr)
    else:
        worker = partial(process_one_frag, args=args)
        with Pool(workers) as pool:
            for msg in pool.map(worker, frag_files):
                print(msg, file=sys.stderr)


if __name__ == "__main__":
    main()
