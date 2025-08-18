# Use a more modern base image
FROM ubuntu:22.04

# Set a working directory
WORKDIR /opt/cd-hit

# Install necessary packages for building cd-hit and for running the parallel scripts
RUN apt-get -y update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    perl \
    openssh-server \
    && rm -rf /var/lib/apt/lists/*

# Copy the entire cd-hit source code into the image
COPY . .

# Compile cd-hit
RUN make

# Compile cd-hit-auxtools
RUN cd /opt/cd-hit/cd-hit-auxtools && make

# Set the system PATH to include all cd-hit executables and scripts
ENV PATH="/opt/cd-hit:/opt/cd-hit/cd-hit-auxtools:/opt/cd-hit/psi-cd-hit:/opt/cd-hit/scripts:${PATH}"

# Expose the SSH port
EXPOSE 22

# Create a directory for the SSH server to run
RUN mkdir /var/run/sshd

# Set a default command to run the SSH server (for worker nodes)
# This can be overridden by the docker-compose file for the master node
CMD ["/usr/sbin/sshd", "-D"]
