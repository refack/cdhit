#!/usr/bin/env python3
import sys
import argparse
from collections import defaultdict
from itertools import combinations
import re

def format_float(f):
    return int(f) if f == int(f) else round(f, 15)

def main():
    parser = argparse.ArgumentParser(
        description="Calculate the sensitivity and specificity of clusters given a benchmark."
    )
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="Input .clstr file (default: stdin)")
    args = parser.parse_args()

    bench_clusters = defaultdict(list)
    cdhit_clusters = defaultdict(list)
    seqs = []
    seq_map = {}

    readin = False
    cluster_idx = 0
    current_cluster = []

    for line in args.input:
        if line.startswith('>'):
            if readin:
                if len(current_cluster) > 1:
                    cdhit_clusters[cluster_idx] = current_cluster
                    cluster_idx += 1
            current_cluster = []
        else:
            readin = True
            match = re.search(r'>([^|]+)\|\|([^.]+)\.\.\.', line)
            if not match:
                match = re.search(r'>([^\.]+)\.\.\.', line)
                if not match:
                    continue
                seq_id, ben_id = match.group(1), None
            else:
                seq_id, ben_id = match.group(1), match.group(2)

            if seq_id not in seq_map:
                seq_map[seq_id] = len(seqs)
                seqs.append(seq_id)

            seq_idx_val = seq_map[seq_id]
            current_cluster.append(seq_idx_val)
            if ben_id:
                bench_clusters[ben_id].append(seq_idx_val)

    if readin and len(current_cluster) > 1:
        cdhit_clusters[cluster_idx] = current_cluster


    cdhit_pairs = set()
    for cluster in cdhit_clusters.values():
        for pair in combinations(sorted(cluster), 2):
            cdhit_pairs.add(pair)

    bench_pairs = set()
    for cluster in bench_clusters.values():
        if len(cluster) > 1:
            for pair in combinations(sorted(cluster), 2):
                bench_pairs.add(pair)

    correct_pairs = cdhit_pairs.intersection(bench_pairs)

    total_bench_pairs = len(bench_pairs)
    total_cdhit_pairs = len(cdhit_pairs)
    total_correct_pairs = len(correct_pairs)

    sensitivity = total_correct_pairs / total_bench_pairs if total_bench_pairs > 0 else 0
    specificity = total_correct_pairs / total_cdhit_pairs if total_cdhit_pairs > 0 else 0

    print(f"Total benchmark pairs\t{total_bench_pairs}")
    print(f"Total cd-hit pairs\t{total_cdhit_pairs}")
    print(f"Total correct pairs\t{total_correct_pairs}")
    print(f"Sensitivity\t{format_float(sensitivity)}")
    print(f"Specificity\t {format_float(specificity)}")

    only_bench = bench_pairs - cdhit_pairs
    print("\n\nPairs in benchmark but not in cd-hit")
    for p1, p2 in sorted(list(only_bench)):
        print(f"{seqs[p1]}\t{seqs[p2]}")

    only_cdhit = cdhit_pairs - bench_pairs
    print("\n\nPairs in cd-hit but not in benchmark")
    for p1, p2 in sorted(list(only_cdhit)):
        print(f"{seqs[p1]}\t{seqs[p2]}")

    print("\n\nPairs in both cd-hit and benchmark")
    for p1, p2 in sorted(list(correct_pairs)):
        print(f"{seqs[p1]}\t{seqs[p2]}")

if __name__ == "__main__":
    main()
