#!/usr/bin/env python3
# =============================================================================
# CD-HIT
# http://cd-hit.org/
#
# program written by
#                                      Weizhong Li
#                                      UCSD, San Diego Supercomputer Center
#                                      La Jolla, CA, 92093
#                                      Email liwz@sdsc.edu
#
# Rewritten in Python by Jules
# =============================================================================

import sys
import re

def get_rep_name(cluster_lines):
    """Extracts the representative sequence name from a list of cluster lines."""
    for line in cluster_lines:
        if line.strip().endswith('*'):
            # Example line: 0	2799aa, >PF04998.6|RPOC2_CHLRE/275-3073... *
            match = re.search(r'>([^\s\.]+)\.\.\.', line)
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
    if len(sys.argv) < 2:
        print("Usage: clstr_merge.py master.clstr slave1.clstr [slave2.clstr ...]", file=sys.stderr)
        sys.exit(1)

    master_clstr_file = sys.argv[1]
    slave_clstr_files = sys.argv[2:]

    slave_file_handles = [open(f, 'r') for f in slave_clstr_files]

    with open(master_clstr_file, 'r') as master_fh:
        for master_cluster in read_clusters(master_fh):
            master_rep_name = get_rep_name(master_cluster)

            # Print the master cluster header and its members
            sys.stdout.write("".join(master_cluster))

            merged_seq_count = len(master_cluster) - 1

            if not master_rep_name:
                continue

            # Process slave files
            for fh in slave_file_handles:
                fh.seek(0) # Rewind the file to scan from the beginning
                for slave_cluster in read_clusters(fh):
                    slave_rep = get_rep_name(slave_cluster)
                    if slave_rep == master_rep_name:
                        # Found the matching cluster in the slave file
                        for line in slave_cluster[1:]:
                            if not line.strip().endswith('*'): # only add non-representatives
                                parts = line.split('\t', 1)
                                if len(parts) == 2:
                                    sys.stdout.write(f"{merged_seq_count}\t{parts[1]}")
                                merged_seq_count += 1
                        break # Move to the next slave file

    for fh in slave_file_handles:
        fh.close()

if __name__ == "__main__":
    main()
