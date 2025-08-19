#!/usr/bin/env python3
import sys
import argparse
import re
from collections import defaultdict

def get_rep_name(cluster_lines):
    """Extracts the representative sequence name from a list of cluster lines."""
    for line in cluster_lines:
        if line.strip().endswith('*'):
            match = re.search(r'>(\S+)\.\.\.', line)
            if match:
                return match.group(1)
    return None

def read_clusters(file_handle):
    """A generator that reads a .clstr file and yields one cluster at a time."""
    cluster_lines = []
    for line in file_handle:
        if line.startswith('>'):
            if cluster_lines:
                yield cluster_lines
            cluster_lines = [line]
        else:
            cluster_lines.append(line)
    if cluster_lines:
        yield cluster_lines

def main():
    parser = argparse.ArgumentParser(
        description="Merge slave .clstr files into a master .clstr file. The order of clusters does not need to be the same."
    )
    parser.add_argument('master_clstr', help="Master .clstr file")
    parser.add_argument('slave_clstrs', nargs='+', help="Slave .clstr files")
    args = parser.parse_args()

    slave_clusters = defaultdict(list)
    for slave_file in args.slave_clstrs:
        with open(slave_file, 'r') as f:
            for cluster in read_clusters(f):
                rep_name = get_rep_name(cluster)
                if rep_name:
                    # store non-representative members
                    for line in cluster:
                        if not line.strip().endswith('*'):
                            slave_clusters[rep_name].append(line)

    with open(args.master_clstr, 'r') as master_fh:
        for master_cluster in read_clusters(master_fh):
            sys.stdout.write("".join(master_cluster))
            master_rep_name = get_rep_name(master_cluster)

            if not master_rep_name:
                continue

            master_members = [l for l in master_cluster if not l.startswith('>')]
            merged_seq_count = len(master_members)

            if master_rep_name in slave_clusters:
                for slave_line in slave_clusters[master_rep_name]:
                    parts = slave_line.split('\t', 1)
                    if len(parts) == 2:
                        sys.stdout.write(f"{merged_seq_count}\t{parts[1]}")
                        merged_seq_count += 1

if __name__ == "__main__":
    main()
