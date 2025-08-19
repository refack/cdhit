#!/usr/bin/env python3
import sys
import argparse
from collections import defaultdict
import re

def format_float(f):
    return int(f) if f == int(f) else f

def main():
    parser = argparse.ArgumentParser(
        description="Calculate the sensitivity and specificity of clusters given a benchmark, using links."
    )
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="Input .clstr file (default: stdin)")
    args = parser.parse_args()

    bench_clusters = defaultdict(list)
    total_cdhit_links = 0
    correct_links = 0

    lines = args.input.readlines()

    # First pass to build bench_clusters from the whole file
    all_seqs = []
    for line in lines:
        if not line.startswith('>'):
            match = re.search(r'>([^|]+)\|\|([^.]+)\.\.\.', line)
            if match:
                seq_id, ben_id = match.group(1), match.group(2)
                all_seqs.append({'id': seq_id, 'ben_id': ben_id})

    for seq in all_seqs:
        bench_clusters[seq['ben_id']].append(seq['id'])

    total_bench_links = sum(len(c) - 1 for c in bench_clusters.values() if len(c) > 1)

    # Second pass to process clusters one by one
    clstr_by_ben = defaultdict(int)
    t_no = 0
    readin = False

    for line in lines:
        if line.startswith('>'):
            if readin and t_no > 1:
                for ben_id in clstr_by_ben:
                    correct_links += clstr_by_ben[ben_id] - 1
                total_cdhit_links += t_no - 1

            t_no = 0
            clstr_by_ben = defaultdict(int)
        else:
            readin = True
            match = re.search(r'>([^|]+)\|\|([^.]+)\.\.\.', line)
            if match:
                ben_id = match.group(2)
                clstr_by_ben[ben_id] += 1
                t_no += 1

    if readin and t_no > 1:
        for ben_id in clstr_by_ben:
            correct_links += clstr_by_ben[ben_id] - 1
        total_cdhit_links += t_no - 1

    sensitivity = correct_links / total_bench_links if total_bench_links > 0 else 0
    specificity = correct_links / total_cdhit_links if total_cdhit_links > 0 else 0

    print(f"Total benchmark links\t{total_bench_links}")
    print(f"Total cd-hit links\t{total_cdhit_links}")
    print(f"Total correct links\t{correct_links}")
    print(f"Sensitivity\t{format_float(sensitivity)}")
    print(f"Specificity\t {format_float(specificity)}")

if __name__ == "__main__":
    main()
