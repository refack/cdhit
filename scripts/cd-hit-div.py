#!/usr/bin/env python3
import sys

def main():
    if len(sys.argv) != 4:
        print("Usage: cd-hit-div.py <input_file> <output_prefix> <num_divisions>", file=sys.stderr)
        sys.exit(1)

    input_file = sys.argv[1]
    output_prefix = sys.argv[2]
    try:
        num_divisions = int(sys.argv[3])
    except ValueError:
        print("Error: num_divisions must be an integer", file=sys.stderr)
        sys.exit(1)

    try:
        output_files = [open(f"{output_prefix}-{i}", "w") for i in range(num_divisions)]
    except IOError as e:
        print(f"Error opening output files: {e}", file=sys.stderr)
        sys.exit(1)

    seq_count = 0
    try:
        with open(input_file, "r") as f_in:
            current_seq_lines = []
            for line in f_in:
                if line.startswith(">"):
                    if current_seq_lines:
                        file_index = seq_count % num_divisions
                        output_files[file_index].writelines(current_seq_lines)
                        seq_count += 1
                    current_seq_lines = [line]
                else:
                    current_seq_lines.append(line)

            if current_seq_lines:
                file_index = seq_count % num_divisions
                output_files[file_index].writelines(current_seq_lines)

    except IOError as e:
        print(f"Error reading input file '{input_file}': {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        for f in output_files:
            f.close()

if __name__ == "__main__":
    main()
