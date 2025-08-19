#!/usr/bin/env python3
import sys
import argparse
import re

def main():
    parser = argparse.ArgumentParser(
        description="Extract representative sequence and size for each cluster."
    )
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="Input .clstr file (default: stdin)")
    args = parser.parse_args()

    cluster_id = ""
    rep_id = ""
    cluster_size = 0

    for line in args.input:
        if line.startswith('>'):
            if cluster_size > 0:
                print(f"{cluster_id}\t{rep_id}\t{cluster_size}")

            match = re.search(r'>Cluster (\d+)', line) or re.search(r'>Clstr (\d+)', line)
            if match:
                cluster_id = match.group(1)
            else:
                cluster_id = ""
            rep_id = ""
            cluster_size = 0
        else:
            cluster_size += 1
            if line.strip().endswith('*'):
                match = re.search(r'>([^\.]+)\.\.\.', line)
                if match:
                    rep_id = match.group(1)

    if cluster_size > 0:
        print(f"{cluster_id}\t{rep_id}\t{cluster_size}")

if __name__ == "__main__":
    main()
