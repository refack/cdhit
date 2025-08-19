#!/usr/bin/env python3
import sys
import argparse

def parse_segs(segs_str):
    """Parses the segments string into a mapping from size to segment index."""
    no2seg_idx = {}
    seg_idx = 0
    for seg in segs_str.split(','):
        if '-' in seg:
            b, e = map(int, seg.split('-'))
            for i in range(b, e + 1):
                no2seg_idx[i] = seg_idx
        else:
            no2seg_idx[int(seg)] = seg_idx
        seg_idx += 1
    return no2seg_idx, seg_idx

def main():
    parser = argparse.ArgumentParser(
        description="Reduce the size of a .clstr file by picking one cluster every 'reduce_rate' clusters for different size segments."
    )
    parser.add_argument('clstr_file', help="Input .clstr file")
    parser.add_argument('segs', help="Comma-separated list of size ranges, e.g., '1,2,3-5,6-10'")
    parser.add_argument('reduce_rate', type=int, help="Reduction rate")
    args = parser.parse_args()

    no2seg_idx, seg_count = parse_segs(args.segs)
    segs_no = [0] * seg_count

    with open(args.clstr_file, 'r') as f:
        cluster_lines = []
        for line in f:
            if line.startswith('>'):
                if cluster_lines:
                    cluster_size = len(cluster_lines) -1
                    if cluster_size in no2seg_idx:
                        this_seg = no2seg_idx[cluster_size]
                        if segs_no[this_seg] % args.reduce_rate == 0:
                            sys.stdout.write("".join(cluster_lines))
                        segs_no[this_seg] += 1
                cluster_lines = [line]
            else:
                cluster_lines.append(line)

        if cluster_lines:
            cluster_size = len(cluster_lines) - 1
            if cluster_size in no2seg_idx:
                this_seg = no2seg_idx[cluster_size]
                if segs_no[this_seg] % args.reduce_rate == 0:
                    sys.stdout.write("".join(cluster_lines))
                segs_no[this_seg] += 1

if __name__ == "__main__":
    main()
