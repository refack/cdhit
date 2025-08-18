# Distributed CD-HIT Execution with Docker Compose

This guide explains how to run `cd-hit` in a distributed fashion across multiple Docker containers. It uses the `cd-hit-para.pl` script, which is part of the standard CD-HIT distribution.

The core challenge in running `cd-hit-para.pl` in a containerized environment is that it assumes a shared filesystem is available to all nodes (master and workers) at the same path. This setup uses Docker Compose and a named volume to create this shared filesystem, allowing the script to work with minimal modification.

## Architecture

The `docker-compose.yml` file defines the following services:

- `master`: The container where you will execute the `cd-hit-para.pl` script. It has SSH client installed to connect to the workers.
- `worker1`, `worker2`: Containers that run an SSH server and wait for commands from the `master`. They perform the actual clustering tasks.
- `shared-data`: A named Docker volume that is mounted to `/data` in all containers. This acts as our shared filesystem.
- `cdhit-net`: A network that allows the containers to communicate with each other by their service names (e.g., `master`, `worker1`).

## Prerequisites

- Docker
- Docker Compose

## Step-by-Step Instructions

### 1. Build the Docker Image

First, build the Docker image for the `cd-hit` environment. All services (master and workers) will use this same image.

Navigate to the `distributed_run` directory and run:

```bash
docker-compose build
```

This command will:
1.  Read the `docker-compose.yml` file.
2.  Find the `build` configuration, which points to the project's root directory (`..`) and the `Docker/Dockerfile`.
3.  Build the image, which includes compiling the `cd-hit` source code from the repository.

### 2. Start the Services

Once the image is built, start all the services in the background:

```bash
docker-compose up -d
```

You can check that the containers are running with `docker-compose ps`. You should see `master`, `worker1`, and `worker2` listed as "up".

### 3. Set Up Passwordless SSH

The `cd-hit-para.pl` script uses SSH to send commands from the `master` to the `worker` nodes. For this to work without manual intervention, you need to set up passwordless SSH from the `master` to each `worker`.

**a. Generate an SSH key on the `master` container:**

```bash
docker-compose exec master ssh-keygen -t rsa -N "" -f /root/.ssh/id_rsa
```
*This command executes `ssh-keygen` inside the `master` container. It creates a new RSA key and saves it to `/root/.ssh/id_rsa` without a passphrase.*

**b. Copy the public key to each worker's `authorized_keys` file:**

For `worker1`:
```bash
docker-compose exec master sh -c "cat /root/.ssh/id_rsa.pub | ssh worker1 'mkdir -p /root/.ssh && cat >> /root/.ssh/authorized_keys'"
```

For `worker2`:
```bash
docker-compose exec master sh -c "cat /root/.ssh/id_rsa.pub | ssh worker2 'mkdir -p /root/.ssh && cat >> /root/.ssh/authorized_keys'"
```
*These commands might prompt you to accept the host key. Just type "yes" and press Enter.*

**c. (Optional) Test the connection:**

You can verify that passwordless SSH is working by trying to connect to a worker:
```bash
docker-compose exec master ssh worker1 "echo 'Hello from worker1'"
```
This should execute the command on `worker1` without asking for a password.

### 4. Run a Distributed CD-HIT Job

Now you are ready to run a distributed clustering job.

**a. Copy data to the shared volume:**

The input data needs to be in the shared volume (`/data`). We have provided a test file `test.fas` and a `hosts` file in this directory. We will copy them into the running `master` container's `/data` directory, which is the shared volume.

```bash
docker cp test.fas $(docker-compose ps -q master):/data/test.fas
docker cp hosts $(docker-compose ps -q master):/data/hosts
```

**b. Execute the `cd-hit-para.pl` script:**

Enter an interactive shell in the `master` container to run the command.

```bash
docker-compose exec master /bin/bash
```

Now, inside the `master` container's shell, run the parallel script. All paths should be relative to the shared `/data` directory.

```bash
cd /data
cd-hit-para.pl -i test.fas -o test.output --B hosts --S 2 -c 0.9
```

**Explanation of the command:**
- `-i test.fas`: The input FASTA file.
- `-o test.output`: The output file for the representative sequences. The cluster file will be `test.output.clstr`.
- `--B hosts`: The file containing the list of worker hostnames.
- `--S 2`: The number of segments to split the data into. This should generally be much larger, but for this small example, 2 is sufficient.
- `-c 0.9`: A standard `cd-hit` option for a 90% sequence identity threshold.

You will see the script print its progress as it sends jobs to `worker1` and `worker2`. When it's finished, you can find the output files (`test.output` and `test.output.clstr`) in the `/data` directory. Since `/data` is a mounted volume, these files will persist even after the containers are stopped.

### 5. Clean Up

Once you are finished, you can stop and remove the containers and network:

```bash
docker-compose down
```

To also remove the shared volume (and all data within it), run:
```bash
docker-compose down -v
```
