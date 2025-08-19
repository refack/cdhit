#!/usr/bin/env python3
import sys
import argparse
import pickle

def get_cluster_size(cluster):
    """Recursively calculates the size of a cluster."""
    size = 0
    # cluster is a list of sequences
    # seq_record: [id, length, has_sub_cluster, sub_cluster, identity]
    for seq in cluster:
        if seq[2]: # has_sub_cluster
            size += get_cluster_size(seq[3])
        else:
            size += 1
    return size

def main():
    parser = argparse.ArgumentParser(
        description="Read a store file, sort the clusters, and save them to a new store file."
    )
    parser.add_argument('input_file', help="Input store file (pickle format)")
    parser.add_argument('output_file', help="Output store file (pickle format)")
    parser.add_argument('sort_by', nargs='?', choices=['no', 'len', 'des'], default='no',
                        help="Sort by: 'no' (number of sequences), 'len' (length of representative), or 'des' (description of representative)")
    args = parser.parse_args()

    with open(args.input_file, 'rb') as f:
        clstr_dict = pickle.load(f)

    clstr_list = list(clstr_dict.values())

    if args.sort_by == 'no':
        # clstr_record: [rep_id, rep_len, has_sub_cluster, sequences, ""]
        clstr_list.sort(key=lambda c: get_cluster_size(c[3]), reverse=True)
    elif args.sort_by == 'len':
        clstr_list.sort(key=lambda c: c[1], reverse=True)
    elif args.sort_by == 'des':
        clstr_list.sort(key=lambda c: c[0])

    with open(args.output_file, 'wb') as f:
        pickle.dump(clstr_list, f)

if __name__ == "__main__":
    main()
