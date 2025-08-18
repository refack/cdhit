#!/usr/bin/env python3
import argparse
from pathlib import Path

def parse_segments(segs_str):
    """Parses the segment string (e.g., '1-3,5') into a list of segment indices."""
    no2seg_idx = {}
    segs = segs_str.split(',')
    for i, seg in enumerate(segs):
        if '-' in seg:
            b, e = map(int, seg.split('-'))
            for j in range(b, e + 1):
                no2seg_idx[j] = i
        else:
            b = int(seg)
            no2seg_idx[b] = i
    return no2seg_idx, len(segs)

def reduce_clusters(input_file, segs_str, reduce_rate):
    """
    Reduces the number of clusters based on a reduction rate per segment.
    """
    no2seg_idx, num_segs = parse_segments(segs_str)
    segs_no = [0] * num_segs

    with open(input_file, 'r') as f:
        this_clstr = ''
        this_no = 0
        readin = False
        for line in f:
            if line.startswith('>'):
                if readin:
                    this_seg = no2seg_idx.get(this_no)
                    if this_seg is not None and (segs_no[this_seg] % reduce_rate) == 0:
                        print(this_clstr, end='')
                    if this_seg is not None:
                        segs_no[this_seg] += 1

                this_no = 0
                this_clstr = line
                readin = False # Will be set to true when we read the first sequence line
            else:
                this_clstr += line
                if not readin:
                    # This logic assumes the number of sequences is what matters.
                    # The original script increments this_no for each sequence line.
                    this_no += 1
                    readin = True

        # Process the last cluster in the file
        if readin:
            this_seg = no2seg_idx.get(this_no)
            if this_seg is not None and (segs_no[this_seg] % reduce_rate) == 0:
                print(this_clstr, end='')

def main():
    parser = argparse.ArgumentParser(description='Reduce CD-HIT cluster file.')
    parser.add_argument('input_file', type=Path, help='Input cluster file (from cd-hit with -d > 0)')
    parser.add_argument('segments', help='Comma-separated list of segment ranges (e.g., "1-10,15,20-25")')
    parser.add_argument('reduce_rate', type=int, help='Reduction rate (integer)')
    args = parser.parse_args()

    if not args.input_file.is_file():
        raise FileNotFoundError(f"Input file not found: {args.input_file}")

    reduce_clusters(args.input_file, args.segments, args.reduce_rate)

if __name__ == '__main__':
    main()

