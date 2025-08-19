#!/usr/bin/env python3
import sys
import argparse
import re

def main():
    parser = argparse.ArgumentParser(
        description="Renumber clusters and sequences in a .clstr file."
    )
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="Input .clstr file (default: stdin)")
    args = parser.parse_args()

    cluster_no = 0
    seq_no = 0
    for line in args.input:
        if line.startswith('>'):
            print(f">Cluster {cluster_no}")
            cluster_no += 1
            seq_no = 0
        else:
            # Replace the first number (sequence index) with the new one
            print(re.sub(r'^\d+', str(seq_no), line), end='')
            seq_no += 1

if __name__ == "__main__":
    main()
