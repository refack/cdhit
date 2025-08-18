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

import argparse
import os
import sys
import subprocess
import json
import shutil

def assign_commands(cd_hit_div_exe, cd_hit_exe, cd_hit_2d_exe, args, arg_pass, indiv):
    commands = []

    # First command is to split the database
    abs_input_path = os.path.abspath(args.i)
    cmd = f"{sys.executable} {cd_hit_div_exe} {abs_input_path} {indiv} {args.S}"
    div_outputs = [f"{indiv}-{i}" for i in range(args.S)]
    commands.append({"cmd": cmd, "status": "wait", "in": [args.i], "out": div_outputs})

    for i in range(args.S):
        idb_orig = f"{indiv}-{i}"
        idb = idb_orig
        idblog = f"{indiv}-{i}.log"

        # compare to previous segs
        for j in range(i):
            jdb = f"{indiv}-{j}-o"
            idbo = f"{indiv}-{i}.vs.{j}"
            cmd = f"{cd_hit_2d_exe} -i {jdb} -i2 {idb} -o {idbo} {arg_pass} >> {idblog}"
            commands.append({"cmd": cmd, "status": "wait", "in": [jdb, idb], "out": [idbo, f"{idbo}.clstr"]})
            idb = idbo

        # self comparing
        idbo = f"{indiv}-{i}-o"
        cmd = f"{cd_hit_exe} -i {idb} -o {idbo} {arg_pass} >> {idblog}"
        in_files = [idb]
        if idb != idb_orig:
            in_files.append(f"{idb}.clstr")
        commands.append({"cmd": cmd, "status": "wait", "in": in_files, "out": [idbo, f"{idbo}.clstr"]})

    return commands

def renumber_clstr(input_file, output_file):
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        cluster_no = 0
        seq_no = 0
        for line in fin:
            if line.startswith('>Cluster'):
                fout.write(f">Cluster {cluster_no}\n")
                cluster_no += 1
                seq_no = 0
            else:
                parts = line.split('\t')
                parts[0] = str(seq_no)
                fout.write('\t'.join(parts))
                seq_no += 1

def main():
    parser = argparse.ArgumentParser(
        description="This script divides a big clustering job into pieces and runs them. This simplified version runs them sequentially.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-i', required=True, help="input filename in fasta format, required")
    parser.add_argument('-o', required=True, help="output filename, required")
    parser.add_argument('--P', default='cd-hit', help='program, "cd-hit" or "cd-hit-est", default "cd-hit"')
    parser.add_argument('--S', type=int, default=64, help="Number of segments to split input DB into, default 64")
    # Ignored arguments for compatibility
    parser.add_argument('--B', help="Ignored")
    parser.add_argument('--L', type=int, default=0, help="Ignored")
    parser.add_argument('--Q', type=int, default=0, help="Ignored")
    parser.add_argument('--T', default='PBS', help='Ignored')
    parser.add_argument('--R', help="Ignored")

    args, unknown = parser.parse_known_args()

    cd_hit_exe = 'cd-hit'
    cd_hit_2d_exe = 'cd-hit-2d'
    if args.P == 'cd-hit-est':
        cd_hit_exe = 'cd-hit-est'
        cd_hit_2d_exe = 'cd-hit-est-2d'

    script_dir = os.path.dirname(os.path.realpath(__file__))
    cd_hit_div_exe = os.path.join(script_dir, 'cd-hit-div.py')
    clstr_merge_exe = os.path.join(script_dir, 'clstr_merge.py')

    arg_pass = ' '.join(unknown)

    work_dir = f"{args.o}.cd-hit-para-tmp"
    indiv = os.path.join(work_dir, os.path.basename(args.i) + ".div")

    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)

    commands = assign_commands(cd_hit_div_exe, cd_hit_exe, cd_hit_2d_exe, args, arg_pass, indiv)

    print("Running commands sequentially...")
    for job in commands:
        print(f"Executing: {job['cmd']}")
        try:
            subprocess.run(job['cmd'], shell=True, check=True, timeout=300)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            print(f"Error executing command: {job['cmd']}\n{e}", file=sys.stderr)
            sys.exit(1)

    print("\nExecution finished. Merging results...")

    out_clstr_merged = f"{args.o}.clstr.merged"
    with open(out_clstr_merged, 'w') as f_out:
        reps_to_cat = []
        for i in range(args.S):
            master_clstr = f"{indiv}-{i}-o.clstr"
            reps_to_cat.append(f"{indiv}-{i}-o")

            slave_clstrs = []
            for j in range(i + 1, args.S):
                tclstr = f"{indiv}-{j}.vs.{i}.clstr"
                if os.path.exists(tclstr):
                    slave_clstrs.append(tclstr)

            if slave_clstrs:
                merge_cmd = [sys.executable, clstr_merge_exe, master_clstr] + slave_clstrs
                result = subprocess.run(merge_cmd, capture_output=True, text=True)
                f_out.write(result.stdout)
            else:
                with open(master_clstr, 'r') as f_master:
                    f_out.write(f_master.read())

    print("Renumbering cluster file...")
    out_clstr = f"{args.o}.clstr"
    renumber_clstr(out_clstr_merged, out_clstr)

    print("Concatenating representative sequences...")
    with open(args.o, 'wb') as f_out:
        for rep_file in reps_to_cat:
            with open(rep_file, 'rb') as f_in:
                shutil.copyfileobj(f_in, f_out)

    print("Done.")

if __name__ == '__main__':
    main()
