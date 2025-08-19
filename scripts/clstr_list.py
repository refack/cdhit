#!/usr/bin/env python3
import sys
import argparse
import pickle
from cdhit.utils import clstr_iter

def main():
    parser = argparse.ArgumentParser(
        description="Parse a .clstr file and store the cluster information in a pickle file. "
                    "If the store file already exists, it merges the new cluster information with the old one."
    )
    parser.add_argument('clstr_file', help="Input .clstr file")
    parser.add_argument('store_file', help="Output store file (pickle format)")
    args = parser.parse_args()

    clusters = {}
    with open(args.clstr_file, 'r') as f:
        for cluster_data in clstr_iter(f):
            rep = None
            for seq in cluster_data:
                if seq['is_rep']:
                    rep = seq
                    break
            if not rep:
                continue

            # seq_record: [id, length, has_sub_cluster, sub_cluster, identity]
            sequences = []
            for seq in cluster_data:
                sequences.append([seq['id'], seq['length'], 0, [], seq['identity']])

            # cluster_record: [rep_id, rep_len, has_sub_cluster, sequences, ""]
            clusters[rep['id']] = [rep['id'], rep['length'], 1, sequences, ""]

    try:
        with open(args.store_file, 'rb') as f:
            old_clstr = pickle.load(f)

        for rep_acc, clstr_data in clusters.items():
            seqs = clstr_data[3]
            for i, seq in enumerate(seqs):
                if seq[0] in old_clstr:
                    clstr_data[3][i][3] = old_clstr[seq[0]][3]
                    clstr_data[3][i][2] = 1
    except FileNotFoundError:
        pass # No old store file, so nothing to merge

    with open(args.store_file, 'wb') as f:
        pickle.dump(clusters, f)

if __name__ == "__main__":
    main()
