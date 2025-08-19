#!/usr/bin/env python3
import sys
import argparse
from cdhit.utils import clstr_iter


def main():
    parser = argparse.ArgumentParser(description="Convert a .clstr file to a tab-separated text file.")
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="Input .clstr file (default: stdin)")
    args = parser.parse_args()

    print("id\tclstr\tclstr_size\tlength\tclstr_rep\tclstr_iden\tclstr_cov")

    for cluster in clstr_iter(args.input):
        cluster.sort(key=lambda x: (x['is_rep'], x['length']), reverse=True)
        rep_len = 0
        for seq in cluster:
            if seq['is_rep']:
                rep_len = seq['length']
                break

        if rep_len == 0: # fallback if no rep found
            rep_len = cluster[0]['length'] if cluster else 1


        for seq in cluster:
            coverage = int(seq['length'] / rep_len * 100) if rep_len else 0
            print(f"{seq['id']}\t"
                  f"{seq['cluster_no']}\t"
                  f"{len(cluster)}\t"
                  f"{seq['length']}\t"
                  f"{int(seq['is_rep'])}\t"
                  f"{seq['identity']}\t"
                  f"{coverage}%")


if __name__ == "__main__":
    main()
