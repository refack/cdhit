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
import sys
import subprocess
import time
import shlex
import shutil
import logging
from pathlib import Path

import dask
from dask.distributed import Client, LocalCluster, as_completed
from dask_jobqueue import PBSCluster, SGECluster


def renumber_clstr(input_file, output_file):
    logging.info(f"Renumbering cluster file {input_file} to {output_file}")
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        cluster_no = 0
        seq_no = 0
        for line in fin:
            if line.startswith('>Cluster'):
                fout.write(f'>Cluster {cluster_no}\n')
                cluster_no += 1
                seq_no = 0
            else:
                parts = line.split('\t')
                parts[0] = str(seq_no)
                fout.write('\t'.join(parts))
                seq_no += 1


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='This script divides a big clustering job into pieces and submits jobs to remote computers over a network to make it parallel.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-i', '--input-file', required=True, help='input filename in fasta format, required')
    parser.add_argument('-o', '--output-directory', required=True, help='output directory, required')
    parser.add_argument('-P', '--program', default='cd-hit-est',
                        help='program, "cd-hit" or "cd-hit-est", default "cd-hit-est"')
    parser.add_argument('-B', '--hosts-file',
                        help='filename of list of hosts, required unless -Q or -L option is supplied')
    parser.add_argument('-L', '--local-cpus', type=int, default=0, help='number of cpus on local computer, default 0')
    parser.add_argument('-S', '--segments', type=int, default=64,
                        help='Number of segments to split input DB into, default 64')
    parser.add_argument('-Q', '--queue-jobs', type=int, default=0,
                        help='number of jobs to submit to queue queuing system, default 0')
    parser.add_argument('-T', '--queue-type', default='PBS',
                        help='type of queuing system, "PBS", "SGE" are supported, default PBS')
    parser.add_argument('--dask-scheduler', help='Address of a Dask scheduler')
    return parser.parse_known_args()


def setup_environment(args):
    """Sets up executables and working directories."""
    logging.info("Setting up environment")
    script_dir = Path(__file__).parent.resolve()
    root_dir = script_dir.parent  # Go up one level from scripts/ to root


    cd_hit_stem = 'cdhit.ps1'
    cd_hit_exe = f'{cd_hit_stem} "" --%'
    cd_hit_2d_exe = f'{cd_hit_stem} 2d --%'
    if args.program == 'cd-hit-est':
        cd_hit_exe = f'{cd_hit_stem} est --%'
        cd_hit_2d_exe = f'{cd_hit_stem} est-2d --%'

    cd_hit_div_exe = script_dir / 'cd-hit-div.py'

    work_dir = Path.cwd() / f'{args.output_directory}.cd-hit-para-tmp'
    indiv = work_dir / Path(args.input_file).stem

    logging.info(f"cd_hit_exe: {cd_hit_exe}")
    logging.info(f"cd_hit_2d_exe: {cd_hit_2d_exe}")
    logging.info(f"cd_hit_div_exe: {cd_hit_div_exe}")
    logging.info(f"work_dir: {work_dir}")
    logging.info(f"indiv: {indiv}")

    return cd_hit_exe, cd_hit_2d_exe, str(cd_hit_div_exe), work_dir, indiv


def setup_dask_cluster(args, work_dir):
    """Sets up and returns a Dask client."""
    logging.info("Setting up Dask cluster")
    if args.dask_scheduler:
        logging.info(f"Connecting to Dask scheduler at {args.dask_scheduler}")
        return Client(args.dask_scheduler)

    if args.local_cpus > 0:
        logging.info(f"Setting up LocalCluster with {args.local_cpus} workers")
        cluster = LocalCluster(n_workers=args.local_cpus, threads_per_worker=1, local_directory=str(work_dir))
        return Client(cluster)

    if args.queue_jobs > 0:
        cluster_class = {'PBS': PBSCluster, 'SGE': SGECluster}.get(args.queue_type)
        if not cluster_class:
            logging.error(f"Unsupported queue type '{args.queue_type}'")
            sys.exit(f"Error: Unsupported queue type '{args.queue_type}'")

        logging.info(f"Setting up {args.queue_type}Cluster with {args.queue_jobs} jobs")
        cluster = cluster_class(
            cores=1,
            memory='2GB',  # Adjust as needed
            n_workers=args.queue_jobs,
            log_directory=str(work_dir / 'dask-logs'),
            local_directory=str(work_dir)
        )
        return Client(cluster)

    if args.hosts_file:
        logging.warning("--hosts-file is not directly supported with Dask in the same way.")
        logging.info("Starting a local cluster. Please start Dask workers on your hosts and connect to the scheduler.")
        cluster = LocalCluster(n_workers=0)  # No local workers unless specified
        logging.info(f'Dask scheduler running at: {cluster.scheduler_address}')
        return Client(cluster)

    return None


def run_command(cmd, cwd, log_file_path):
    """A generic function to run a command and log its output."""
    cmd_str = ' '.join(map(str, cmd))
    logging.info(f"Executing command: {cmd_str}")
    logging.info(f"  CWD: {cwd}")
    logging.info(f"  Log file: {log_file_path}")

    try:
        with open(log_file_path, 'a') as log_file:
            log_file.write(f"Executing command: {cmd_str}\n")
            log_file.write(f"CWD: {cwd}\n")
            result = subprocess.run(cmd, cwd=cwd, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            log_file.write("--- STDOUT ---\n")
            log_file.write(result.stdout)
            log_file.write("--- STDERR ---\n")
            log_file.write(result.stderr)
            logging.info(f"Command finished successfully: {cmd_str}")
            if result.stdout:
                logging.debug(f"STDOUT: {result.stdout.strip()}")
            if result.stderr:
                logging.debug(f"STDERR: {result.stderr.strip()}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed with exit code {e.returncode}: {cmd_str}")
        logging.error(f"  See log for details: {log_file_path}")
        raise
    except Exception as e:
        logging.error(f"An unexpected error occurred while running command: {cmd_str}")
        logging.error(f"  Exception: {e}")
        logging.error(f"  See log for details: {log_file_path}")
        raise


def run_cdhit_task(cd_hit_exe, arg_pass, work_dir, input_file, out_prefix):
    """Dask task for running cd-hit on a single segment."""
    out_reps = f"{out_prefix}"
    out_clstr = f"{out_prefix}.clstr"
    log_file = f"{out_prefix}.log"

    if Path(out_reps).exists() and Path(out_clstr).exists():
        logging.info(f"Outputs for {input_file} already exist, skipping.")
        return out_reps, out_clstr

    cmd = [cd_hit_exe, '-i', input_file, '-o', out_reps] + shlex.split(arg_pass)
    run_command(cmd, cwd=work_dir.parent, log_file_path=log_file)
    return out_reps, out_clstr


def run_cdhit_2d_task(cd_hit_2d_exe, arg_pass, work_dir, db1_files, db2_files, out_prefix):
    """Dask task for running cd-hit-2d to combine two results."""
    db1_reps, _ = db1_files
    db2_reps, _ = db2_files

    out_reps = f"{out_prefix}"
    out_clstr = f"{out_prefix}.clstr"
    log_file = f"{out_prefix}.log"

    if Path(out_reps).exists() and Path(out_clstr).exists():
        logging.info(f"Outputs for 2D job {out_prefix} already exist, skipping.")
        return out_reps, out_clstr

    cmd = [cd_hit_2d_exe, '-i', db1_reps, '-i2', db2_reps, '-o', out_reps] + shlex.split(arg_pass)
    run_command(cmd, cwd=work_dir.parent, log_file_path=log_file)
    return out_reps, out_clstr


def execute_tree_reduction(client, cd_hit_exe, cd_hit_2d_exe, cd_hit_div_exe, args, arg_pass, work_dir, indiv):
    """Builds and executes the clustering computation graph using Dask."""
    logging.info("Starting tree reduction execution.")

    # 1. Split database
    logging.info("Step 1: Splitting database")
    div_cmd = [sys.executable, cd_hit_div_exe, str(Path(args.input_file).resolve()), str(indiv), str(args.segments)]
    div_log = work_dir / "cd-hit-div.log"
    segment_files = [work_dir / f"{indiv.name}-{i}" for i in range(args.segments)]
    if not all(f.exists() for f in segment_files):
        run_command(div_cmd, cwd=work_dir, log_file_path=div_log)
    else:
        logging.info("Segment files already exist, skipping database split.")

    # 2. Level 0 - Self-clustering tasks (map)
    logging.info("Step 2: Building level 0 (self-clustering) tasks")
    tasks = []
    for i in range(args.segments):
        input_file = segment_files[i]
        out_prefix = work_dir / f"level-0-seg-{i}"
        task = dask.delayed(run_cdhit_task)(cd_hit_exe, arg_pass, work_dir, input_file, out_prefix, dask_key_name=f"level-0-task-{i}")
        tasks.append(task)

    # 3. Reduction levels
    logging.info("Step 3: Building reduction tree")
    level = 0
    while len(tasks) > 1:
        logging.info(f"Building level {level + 1}: {len(tasks)} tasks -> {len(tasks) // 2} tasks")
        new_tasks = []
        
        # Handle odd number of tasks
        odd_task_out = None
        if len(tasks) % 2 != 0:
            odd_task_out = tasks.pop()

        for i in range(0, len(tasks), 2):
            db1_files_task = tasks[i]
            db2_files_task = tasks[i+1]
            out_prefix = work_dir / f"level-{level + 1}-pair-{i // 2}"
            new_task = dask.delayed(run_cdhit_2d_task)(cd_hit_2d_exe, arg_pass, work_dir, db1_files_task, db2_files_task, out_prefix, dask_key_name=f"level-{level+1}-task-{i//2}")
            new_tasks.append(new_task)
        
        if odd_task_out is not None:
            new_tasks.append(odd_task_out) # Carry over the odd one
        
        tasks = new_tasks
        level += 1

    final_task = tasks[0]

    # 4. Execute graph
    logging.info("Step 4: Computing final result")
    # final_task.visualize(filename='dask-graph.png') # Uncomment for debugging
    final_reps_path, final_clstr_path = final_task.compute()
    logging.info("Tree reduction computation finished.")

    return final_reps_path, final_clstr_path


def merge_results(args, final_reps_path, final_clstr_path):
    """Copies final results to the output directory."""
    logging.info('\nExecution finished. Finalizing results...')

    # Copy final representative sequences file
    logging.info(f"Copying final representative sequences to {args.output_directory}")
    shutil.copy(final_reps_path, args.output_directory)

    # Renumber final cluster file
    out_clstr = f'{args.output_directory}.clstr'
    logging.info(f"Renumbering final cluster file to {out_clstr}")
    renumber_clstr(final_clstr_path, out_clstr)


def main():
    """Main function to run the parallel clustering."""
    args, unknown = parse_arguments()
    arg_pass = ' '.join(unknown)

    log_file = "cd-hit-para.log"
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        handlers=[
                            logging.FileHandler(log_file),
                            logging.StreamHandler(sys.stdout)
                        ])
    
    # Suppress overly verbose logs from Dask
    logging.getLogger('distributed.worker').setLevel(logging.WARNING)
    logging.getLogger('distributed.core').setLevel(logging.WARNING)
    logging.getLogger('distributed.client').setLevel(logging.WARNING)


    logging.info(f"Starting cd-hit-para.py with arguments: {args}")
    logging.info(f"Unknown arguments passed to cd-hit: {arg_pass}")

    cd_hit_exe, cd_hit_2d_exe, cd_hit_div_exe, work_dir, indiv = setup_environment(args)

    # Don't remove work_dir if it exists, to allow for restarts
    if not work_dir.exists():
        logging.info(f"Creating working directory: {work_dir}")
        work_dir.mkdir(parents=True, exist_ok=True)
    else:
        logging.info(f"Working directory {work_dir} already exists. Will reuse existing files.")

    client = setup_dask_cluster(args, work_dir)
    if not client:
        logging.error("Could not set up a Dask cluster. Please specify an execution environment (--local-cpus, --queue-jobs, etc.)")
        sys.exit("Error: Could not set up a Dask cluster.")

    logging.info(f'Dask dashboard link: {client.dashboard_link}')

    try:
        final_reps, final_clstr = execute_tree_reduction(client, cd_hit_exe, cd_hit_2d_exe, cd_hit_div_exe, args, arg_pass, work_dir, indiv)
        merge_results(args, final_reps, final_clstr)
        logging.info('Done.')
    except Exception as e:
        logging.critical(f"A critical error occurred during execution: {e}")
        sys.exit("Exiting due to critical error.")
    finally:
        logging.info("Closing Dask client")
        client.close()


if __name__ == '__main__':
    main()
