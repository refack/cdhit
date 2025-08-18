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

    # First command is to split the database
    abs_input_path = os.path.abspath(args.i)
    cmd = [sys.executable, cd_hit_div_exe, abs_input_path, indiv, str(args.S)]
    div_outputs = [f"{indiv}-{i}" for i in range(args.S)]
    commands.append({"cmd": cmd, "status": "wait", "in": [args.i], "out": div_outputs, "log_idx": -1})

    for i in range(args.S):
        idb_orig = f"{indiv}-{i}"
        idb = idb_orig

        # compare to previous segs
        for j in range(i):
            jdb = f"{indiv}-{j}-o"
            idbo = f"{indiv}-{i}.vs.{j}"
            cmd = [cd_hit_2d_exe, "-i", jdb, "-i2", idb, "-o", idbo] + shlex.split(arg_pass)
            commands.append({"cmd": cmd, "status": "wait", "in": [jdb, idb], "out": [idbo, f"{idbo}.clstr"], "log_idx": i})
            idb = idbo

        # self comparing
        idbo = f"{indiv}-{i}-o"
        cmd = [cd_hit_exe, "-i", idb, "-o", idbo] + shlex.split(arg_pass)
        in_files = [idb]
        if idb != idb_orig:
            in_files.append(f"{idb}.clstr")
        commands.append({"cmd": cmd, "status": "wait", "in": in_files, "out": [idbo, f"{idbo}.clstr"], "log_idx": i})

    return commands

def write_restart(restart_file, commands):
    # Convert command lists to strings for JSON serialization
    for cmd in commands:
        if isinstance(cmd['cmd'], list):
            cmd['cmd_str'] = ' '.join(cmd['cmd'])
        else:
            cmd['cmd_str'] = cmd['cmd']

    with open(restart_file, 'w') as f:
        json.dump([{'cmd_str': c['cmd_str'], 'status': c['status'], 'in': c['in'], 'out': c['out'], 'log_idx': c['log_idx']} for c in commands], f, indent=2)

def read_restart(restart_in):
    with open(restart_in) as f:
        data = json.load(f)
    # Convert command strings back to lists
    for item in data:
        item['cmd'] = shlex.split(item['cmd_str'])
    return data


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
        time.sleep(1) # Reduced sleep time
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
        description="This script divides a big clustering job into pieces and submits jobs to remote computers over a network to make it parallel.",
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

    script_dir = os.path.dirname(os.path.realpath(__file__))
    cd_hit_div_exe = os.path.join(script_dir, 'cd-hit-div.py')
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
    running_jobs = {} # host_idx -> {proc: Popen, job_idx: int}

    if args.R and os.path.exists(args.R):
        print(f"Restarting from {args.R}")
        commands = read_restart(args.R)
    else:
        if os.path.exists(work_dir):
            shutil.rmtree(work_dir)
        os.makedirs(work_dir)
        time.sleep(1)
        commands = assign_commands(cd_hit_div_exe, cd_hit_exe, cd_hit_2d_exe, args, arg_pass, indiv)
        #write_restart(restart_file, commands)

    if commands[0]['status'] == 'wait':
        print(f"Executing: {' '.join(commands[0]['cmd'])}")
        subprocess.run(commands[0]['cmd'], check=True)
        commands[0]['status'] = 'done'
        #write_restart(restart_file, commands)
        for out_file in commands[0]['out']:
            stable_files.add(out_file)

    sleep_time = 1
    while True:
        # Check for finished jobs
        for host_idx, job_info in list(running_jobs.items()):
            proc = job_info['proc']
            job_idx = job_info['job_idx']
            if proc.poll() is not None:
                if proc.returncode == 0:
                    commands[job_idx]['status'] = 'done'
                    print(f"\nJob {job_idx} completed successfully.")
                else:
                    commands[job_idx]['status'] = 'failed'
                    print(f"\nError: Job {job_idx} failed with return code {proc.returncode}. Command: {' '.join(commands[job_idx]['cmd'])}")
                del running_jobs[host_idx]

        all_done = all(job['status'] in ['done', 'failed'] for job in commands)
        if any(job['status'] == 'failed' for job in commands):
            print("One or more jobs failed. Aborting.")
            break
        if all_done:
            print("\nAll jobs completed.")
            break

        # Try to send new jobs
        available_hosts = [i for i in range(host_no) if i not in running_jobs]

        if available_hosts:
            for i, job in enumerate(commands):
                if job['status'] != 'wait':
                    continue

                inputs_ready = all(os.path.exists(f) for f in job.get('in', []))
                if not inputs_ready:
                    continue

                for in_file in job.get('in', []):
                    wait_stable_file(in_file, stable_files)

                host_idx = available_hosts.pop(0)
                host = hosts[host_idx]

                job_cmd_str = ' '.join(job['cmd'])
                print(f"\nExecuting job {i} on {host}: {job_cmd_str}")

                log_idx = job.get('log_idx', -1)
                log_file = None
                if log_idx != -1:
                    log_file_path = os.path.join(work_dir, f"log.{log_idx}")
                    log_file = open(log_file_path, "a")

                proc = None
                if args.L > 0: # Local execution
                    proc = subprocess.Popen(job['cmd'], stdout=log_file or subprocess.DEVNULL, stderr=log_file or subprocess.STDOUT)
                elif args.Q > 0: # Queue submission
                    tsh_file = os.path.join(work_dir, f"job.{i}.sh")
                    with open(tsh_file, 'w') as f:
                        f.write("#!/bin/sh\n")
                        f.write(f"{job_cmd_str}\n")
                    os.chmod(tsh_file, 0o755)

                    qsub_cmd = ""
                    if args.T == 'PBS':
                        qsub_cmd = f"qsub -N para-{i} -o {log_file_path or '/dev/null'} -e {log_file_path or '/dev/null'} {tsh_file}"
                    elif args.T == 'SGE':
                        qsub_cmd = f"qsub -N para-{i} -S /bin/bash -V -cwd -o {log_file_path or '/dev/null'} -e {log_file_path or '/dev/null'} {tsh_file}"

                    if qsub_cmd:
                        proc = subprocess.Popen(shlex.split(qsub_cmd))
                elif args.B: # SSH execution
                    ssh_cmd = f"ssh -xq {host} 'cd {pwd}; {job_cmd_str}'"
                    proc = subprocess.Popen(ssh_cmd, shell=True, stdout=log_file or subprocess.DEVNULL, stderr=log_file or subprocess.STDOUT)

                if proc:
                    running_jobs[host_idx] = {'proc': proc, 'job_idx': i}
                    job['status'] = 'run'

                if log_file:
                    log_file.close()

                if not available_hosts:
                    break

        statuses = [job['status'] for job in commands]
        print(f"\rJob statuses: { {s: statuses.count(s) for s in set(statuses)} }", end="", flush=True)
        time.sleep(1)

    if any(job['status'] == 'failed' for job in commands):
        sys.exit("Exiting due to job failure.")

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
