#!/usr/bin/env python3
"""Recover marker dosages from RTIGER segments on the same coordinate system."""
import argparse
import bisect
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def read_table(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def chromosome(value):
    value = value.lower().removeprefix("chr")
    return "chr" + str(int(value)) if value.isdigit() else value


def recover(markers, calls, segments):
    marker_ids = [m["marker_id"] for m in markers]
    sample_ids = [r["sample_id"] for r in calls]
    if not markers or not calls or len(set(marker_ids)) != len(marker_ids) or len(set(sample_ids)) != len(sample_ids):
        raise ValueError("Marker and sample identifiers must be nonempty and unique")
    groups = defaultdict(list)
    state_to_dosage = {"mat": -1, "het": 0, "pat": 1}
    for row in segments:
        start, end = int(row["first_marker"]), int(row["last_marker"])
        if start < 1 or end < start:
            raise ValueError("Invalid segment coordinates")
        groups[(row["sample"], chromosome(row["chr"]))].append(
            (start, end, state_to_dosage[row["state"]]))
    starts = {}
    for key, spans in groups.items():
        spans.sort()
        if any(left[1] >= right[0] for left, right in zip(spans, spans[1:])):
            raise ValueError(f"Overlapping segments: {key}")
        starts[key] = [s[0] for s in spans]
    result = []
    for row in calls:
        values = [row["sample_id"]]
        for marker in markers:
            key = (row["sample_id"], chromosome(marker["chr"]))
            if key not in groups:
                raise ValueError(f"Missing segments: {key}")
            spans, pos = groups[key], int(marker["pos"])
            if pos < 1:
                raise ValueError("Marker positions are one-based")
            index = bisect.bisect_right(starts[key], pos) - 1
            if index < 0:
                dosage = spans[0][2]
            elif pos <= spans[index][1] or index == len(spans) - 1:
                dosage = spans[index][2]
            else:
                midpoint = (spans[index][1] + spans[index + 1][0]) // 2
                dosage = spans[index if pos <= midpoint else index + 1][2]
            observed = row[marker["marker_id"]]
            if observed not in ("", "NA", ".") and float(observed) != dosage:
                raise ValueError(f"Segment disagrees with observed call: {key}, {marker['marker_id']}")
            values.append(dosage)
        result.append(values)
    return ["sample_id"] + marker_ids, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--markers", required=True, help="TSV: marker_id, chr, pos")
    parser.add_argument("--genotypes", required=True, help="TSV: sample_id followed by marker dosage columns")
    parser.add_argument("--segments", required=True, help="RTIGER raw_viterbi_segments.tsv")
    parser.add_argument("--output", required=True, help="Complete genotype TSV")
    args = parser.parse_args()
    _, markers = read_table(args.markers)
    columns, calls = read_table(args.genotypes)
    _, segments = read_table(args.segments)
    if columns != ["sample_id"] + [m["marker_id"] for m in markers]:
        raise ValueError("Genotype columns must match the marker table in order")
    header, rows = recover(markers, calls, segments)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


if __name__ == "__main__":
    main()
