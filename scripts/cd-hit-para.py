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
import time
import json
import shlex
import shutil

def assign_commands(cd_hit_div_exe, cd_hit_exe, cd_hit_2d_exe, args, arg_pass, indiv):
    commands = []

    cmd = f"{cd_hit_div_exe} -i {args.i} -o {indiv} -div {args.S}"
    outputs = [f"{indiv}-{i}" for i in range(args.S)]
    for out in outputs:
        outputs.append(f"{out}.clstr")
    commands.append({"cmd": cmd, "status": "wait", "in": [args.i], "out": outputs})

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
        if idb != idb_orig: # it was an output of a previous step
            in_files.append(f"{idb}.clstr")
        commands.append({"cmd": cmd, "status": "wait", "in": in_files, "out": [idbo, f"{idbo}.clstr"]})

    return commands

def write_restart(restart_file, commands):
    with open(restart_file, 'w') as f:
        json.dump(commands, f, indent=2)

def read_restart(restart_in):
    with open(restart_in) as f:
        return json.load(f)

def wait_stable_file(f, stable_files):
    if f in stable_files:
        return
    if not os.path.exists(f):
        return

    if os.path.exists(f"{f}.done"):
        stable_files.add(f)
        return

    size0 = os.path.getsize(f)
    while True:
        time.sleep(10)
        size1 = os.path.getsize(f)
        if size0 == size1:
            stable_files.add(f)
            break
        else:
            size0 = size1

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
        description="This script divides a big clustering job into pieces and submits jobs to remote computers over a network to make it parallel. After all the jobs finished, the script merges the clustering results as if you just run a single cd-hit or cd-hit-est.\n\nYou can also use it to divide big jobs on a single computer if your computer does not have enough RAM (with -L option).",
        epilog="More cd-hit/cd-hit-est options can be specified in the command line.\n\nQuestions, bugs, contact Weizhong Li at liwz@sdsc.edu",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-i', required=True, help="input filename in fasta format, required")
    parser.add_argument('-o', required=True, help="output filename, required")
    parser.add_argument('--P', default='cd-hit', help='program, "cd-hit" or "cd-hit-est", default "cd-hit"')
    parser.add_argument('--B', help="filename of list of hosts, required unless -Q or -L option is supplied")
    parser.add_argument('--L', type=int, default=0, help="number of cpus on local computer, default 0")
    parser.add_argument('--S', type=int, default=64, help="Number of segments to split input DB into, default 64")
    parser.add_argument('--Q', type=int, default=0, help="number of jobs to submit to queue queuing system, default 0")
    parser.add_argument('--T', default='PBS', help='type of queuing system, "PBS", "SGE" are supported, default PBS')
    parser.add_argument('--R', help="restart file, used after a crash of run")

    args, unknown = parser.parse_known_args()

    cd_hit_exe = 'cd-hit'
    cd_hit_2d_exe = 'cd-hit-2d'
    if args.P == 'cd-hit-est':
        cd_hit_exe = 'cd-hit-est'
        cd_hit_2d_exe = 'cd-hit-est-2d'

    cd_hit_div_exe = 'cd-hit-div'
    script_dir = os.path.dirname(os.path.realpath(__file__))
    clstr_merge_exe = os.path.join(script_dir, 'clstr_merge.py')

    arg_pass = ' '.join(unknown)

    pwd = os.getcwd()
    work_dir = f"{args.o}.cd-hit-para-tmp"
    restart_file = f"{args.o}.restart"
    indiv = os.path.join(work_dir, os.path.basename(args.i) + ".div")

    commands = []
    stable_files = set()

    hosts = []
    if args.B:
        with open(args.B) as f:
            hosts = [line.strip() for line in f if line.strip()]
    elif args.Q > 0:
        hosts = [f"queue_host.{i}" for i in range(args.Q)]
    elif args.L > 0:
        hosts = [f"localhost.{i}" for i in range(args.L)]

    if not hosts:
        sys.exit("Error: no hosts specified. Use --B, --Q, or --L.")

    host_no = len(hosts)
    running_jobs = {} # host_idx -> process object

    if args.R and os.path.exists(args.R):
        print(f"Restarting from {args.R}")
        commands = read_restart(args.R)
    else:
        if os.path.exists(work_dir):
            shutil.rmtree(work_dir)
        os.makedirs(work_dir)
        commands = assign_commands(cd_hit_div_exe, cd_hit_exe, cd_hit_2d_exe, args, arg_pass, indiv)
        write_restart(restart_file, commands)

    # First command (cd-hit-div) is special, run it synchronously
    if commands[0]['status'] == 'wait':
        print(f"Executing: {commands[0]['cmd']}")
        subprocess.run(commands[0]['cmd'], shell=True, check=True)
        commands[0]['status'] = 'done'
        write_restart(restart_file, commands)
        for i in range(args.S):
            stable_files.add(f"{indiv}-{i}")

    sleep_time = 1
    while True:
        # Check for finished jobs
        for host_idx, proc in list(running_jobs.items()):
            if proc.poll() is not None: # Process has terminated
                del running_jobs[host_idx]

        status_change = False
        all_done = True
        for job in commands:
            if job['status'] == 'run':
                all_done = False
            elif job['status'] == 'done':
                continue
            else: # wait
                all_done = False
                # Check if all output files for this job exist and are non-empty
                outputs_exist = True
                if job.get('out'):
                    for out_file in job['out']:
                        if not os.path.exists(out_file): # or os.path.getsize(out_file) == 0:
                            outputs_exist = False
                            break

                if outputs_exist and job.get('out'):
                    job['status'] = 'done'
                    status_change = True
                    print(f"\nJob completed based on file check: {job['cmd']}")

        if status_change:
            write_restart(restart_file, commands)

        if all_done:
            print("\nAll jobs completed.")
            break

        job_sent = False
        available_hosts = [i for i in range(host_no) if i not in running_jobs]

        if available_hosts:
            for job in commands:
                if job['status'] != 'wait':
                    continue

                inputs_ready = True
                for in_file in job.get('in', []):
                    if not os.path.exists(in_file):
                        inputs_ready = False
                        break
                if not inputs_ready:
                    continue

                for in_file in job.get('in', []):
                    wait_stable_file(in_file, stable_files)

                host_idx = available_hosts.pop(0)
                host = hosts[host_idx]

                tsh_file = os.path.join(work_dir, f"{os.path.basename(args.o)}.{os.getpid()}.{host_idx}.sh")

                with open(tsh_file, 'w') as f:
                    f.write("#!/bin/sh\n")
                    f.write(f"set -e\n")
                    f.write(f"{job['cmd']}\n")
                    for out_file in job.get('out', []):
                         f.write(f"touch {out_file}.done\n")

                os.chmod(tsh_file, 0o755)

                job['status'] = 'run'
                print(f"\nExecuting on {host}: {job['cmd']}")

                if args.L > 0: # Local execution
                    proc = subprocess.Popen(['/bin/sh', tsh_file], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
                    running_jobs[host_idx] = proc
                elif args.Q > 0: # Queue submission
                    qsub_cmd = ""
                    if args.T == 'PBS':
                        qsub_cmd = f"qsub -N para-{host_idx} -o {host}.log -e {host}.err {tsh_file}"
                    elif args.T == 'SGE':
                        qsub_cmd = f"qsub -N para-{host_idx} -S /bin/bash -V -cwd {tsh_file}"

                    if qsub_cmd:
                        proc = subprocess.Popen(shlex.split(qsub_cmd))
                else: # SSH execution
                    ssh_cmd = f"ssh -xq {host} 'cd {pwd}; /bin/sh {tsh_file}'"
                    proc = subprocess.Popen(ssh_cmd, shell=True)

                job_sent = True
                if not available_hosts:
                    break

        if not job_sent:
            print(".", end="", flush=True)
            time.sleep(sleep_time)
            if sleep_time < 30:
                sleep_time += 1

    print("\nExecution finished. Merging results...")

    out_clstr_merged = f"{args.o}.clstr.merged"
    with open(out_clstr_merged, 'w') as f_out:
        reps_to_cat = []
        for i in range(args.S):
            master_clstr = f"{indiv}-{i}-o.clstr"
            reps_to_cat.append(f"{indiv}-{i}-o")

            slave_clstrs = []
            for j in range(i + 1, args.S):
                tclstr = f"{indiv}-{i}.vs.{j}.clstr"
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

    # print(f"Cleaning up temporary directory {work_dir}...")
    # shutil.rmtree(work_dir)
    print("Done.")

if __name__ == '__main__':
    main()
