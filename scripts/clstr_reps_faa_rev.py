#!/usr/bin/env python3
import sys
import argparse
import re
from cdhit.utils import fasta_iter

def main():
    parser = argparse.ArgumentParser(
        description="Filter a FASTA file based on a .clstr file, keeping a certain number of sequences per cluster."
    )
    parser.add_argument('clstr_file', help="Input .clstr file")
    parser.add_argument('fasta_file', help="Input FASTA file")
    parser.add_argument('cutoff', type=int, help="Number of sequences to keep per cluster")
    args = parser.parse_args()

    if args.cutoff <= 0:
        sys.exit("Error: cutoff must be greater than 0")

    ids_to_skip = set()
    with open(args.clstr_file, 'r') as f:
        cluster_members = []
        for line in f:
            if line.startswith('>'):
                if cluster_members:
                    if len(cluster_members) > args.cutoff:
                        for i in range(args.cutoff, len(cluster_members)):
                            ids_to_skip.add(cluster_members[i])
                cluster_members = []
            else:
                match = re.search(r'>(\S+?)\.\.\.', line)
                if match:
                    seq_id = match.group(1)
                    if line.strip().endswith('*'):
                        cluster_members.insert(0, seq_id)
                    else:
                        cluster_members.append(seq_id)

        if cluster_members: # process the last cluster
            if len(cluster_members) > args.cutoff:
                for i in range(args.cutoff, len(cluster_members)):
                    ids_to_skip.add(cluster_members[i])

    fasta_file_handle = sys.stdin if args.fasta_file == '-' else open(args.fasta_file, 'r')

    flag = 0
    for line in fasta_file_handle:
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            header = line[1:]
            seq_id_parts = header.split()
            if seq_id_parts:
                seq_id = seq_id_parts[0]
                if seq_id in ids_to_skip:
                    flag = 0
                else:
                    flag = 1
            else:
                flag = 1 # Keep empty headers

        if flag:
            print(line)

    if args.fasta_file != '-':
        fasta_file_handle.close()

if __name__ == "__main__":
    main()
