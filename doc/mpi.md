# Tips and Tricks for using MPI with IBAMR

## Overview
IBAMR is parallelized with MPI. MPI, the message-passing interface, is the most
popular network communication standard for large-scale (dozens, hundreds, or
more processors) parallelization of scientific programs at the current time.
Since MPI is an industry standard there are many high-quality
implementations. The two best-known implementations are OpenMPI and MPICH.

## Getting Started

### Setting up MPI
Whenever possible, use the copy of MPI provided by the system
administrators. This is usually available through the module system.

MPI is a low-level library in the sense that its performance is heavily
dependent on specialized hardware. In particular, to utilize common hardware on
clusters (such as infiniband), the MPI library must be configured in a
machine-specific way that is not trivial to do. Compiling and using your own MPI
library instead of the one provided by the administrators can result in large
performance problems (e.g., code can run ten times slower).

That being said, it may occasionally be necessary to compile MPI on your own for
some reason or another. Since IBAMR depends on PETSc, we recommend using PETSc
to build MPI with the `--download-mpich` option.

IBAMR uses MPI to create multiple copies of the same program (the 'same
instruction, multiple data' paradigm). Each copy runs in a unique process - the
only way for processes to communicate is over the MPI network.

IBAMR does *not* parallelize over multiple threads. Additionally, IBAMR
applications do not benefit much from hyperthreading (i.e., running multiple
processes on a single core). For the remainder of this document we pretend
hyperthreading does not exist since it is not useful to our applications.

### Running programs with MPI
IBAMR applications are always run with a fixed number of processes. They are
started by running, e.g.,
```
mpiexec -n 16 ./main3d ./input3d
```
to run `main3d` with input file `input3d` with `16` processors.

### Determining the number of processors to use
The performance of many parts of IBAMR (such as the fluid solve) is limited by
the rate at which data can be moved from main memory (RAM) into the CPU. The
speed of this transfer can be measured with the `streams` benchmark. One can run
the `streams` benchmark with PETSc by running
```
make streams
```
in PETSc's root directory. The output looks something like this:
```
6  45603.4653   Rate (MB/s) 3.81156
7  42492.4568   Rate (MB/s) 3.55154
8  48441.6695   Rate (MB/s) 4.04878
```
indicating that we can stream 48.4 GB/s with 8 processors and get a factor of 4
speedup (versus using 1 processor), measured by the time it takes to stream
data. This is typical: performance usually reaches a maximum when only a
fraction of a machine's processors are used.

## Improving performance

### Mapping MPI processes to particular physical resources
Performance is dependent on the way processors are mapped both to *cores* (the
actual physical CPU) and *sockets*. Modern machines typically have multiple
cores. Larger machines, such as workstations and servers, have multiple sockets,
which correspond to groups of cores. My workstation has 28 physical processors
which are divided among two sockets. To get this information on a Linux machine
try looking at the special file `/proc/cpuinfo` (here `physical id` refers to
the socket number).

MPI refers to the process of assigning processes to physical resources as
*mapping*. Mapping is usually done at the core, socket, and node level.

For example, if we run
```
$ mpiexec -np 8 --map-by core ./streams
8  24723.6086   Rate (MB/s) 2.06641
```
i.e., we instruct the machine to map to the first 8 cores instead of mapping
evenly across the two sockets in this machine (the default). This drops to
roughly half of the peak speed because we have maximized bandwidth on a single
socket instead of maximizing bandwidth across two sockets.

Note that we get essentially identical results if we specify a number of
processes equal to the number of physical processors on the machine:
```
$ mpiexec -np 28 --map-by core ./streams
28  48480.7239   Rate (MB/s) 4.05205
$ mpiexec -np 28 --map-by socket ./streams
28  48411.4946   Rate (MB/s) 4.04626
```

The exact binding behavior depends on the version of MPI being used. Here, I
have `openmpi` on my workstation, which defaults to binding processes to
sockets - that is why I get good performance with the default.

One can avoid doing some arithmetic by immediately assigning processes to
sockets:
```
$ mpiexec -map-by ppr:4:socket ./streams
8  48606.2252   Rate (MB/s) 4.06254
```
More generally, the arguments `-map-by ppr:N:<object>` creates `N` times the
number of `<object>`s copies of the application on each node.

### Binding processes to particular physical resources
In MPI terminology, *mapping* refers to the initial allocation of resources
(starting things in certain places) while *binding* refers to adding constraints
that require processes to run in specific places.

### Displaying some additional information about MPI
Try using the `-display-map` option to see where different processes ended up:
```
$ mpiexec --display-map --map-by ppr:4:socket ./streams
 Data for JOB [22109,1] offset 0

 ========================   JOB MAP   ========================

 Data for node: auricle Num slots: 28   Max slots: 0    Num procs: 8
        Process OMPI jobid: [22109,1] App: 0 Process rank: 0
        Process OMPI jobid: [22109,1] App: 0 Process rank: 1
        Process OMPI jobid: [22109,1] App: 0 Process rank: 2
        Process OMPI jobid: [22109,1] App: 0 Process rank: 3
        Process OMPI jobid: [22109,1] App: 0 Process rank: 4
        Process OMPI jobid: [22109,1] App: 0 Process rank: 5
        Process OMPI jobid: [22109,1] App: 0 Process rank: 6
        Process OMPI jobid: [22109,1] App: 0 Process rank: 7

 =============================================================
8  48480.5008   Rate (MB/s) 4.05203
```

Similarly, `-report-bindings` creates a nice plot of where processes end up:
```
[drwells@auricle streams]$ mpiexec -report-bindings --bind-to core --map-by ppr:4:socket ./streams
[auricle:14958] MCW rank 4 bound to socket 1[core 14[hwt 0-1]]: [../../../../../../../../../../../../../..][BB/../../../../../../../../../../../../..]
[auricle:14958] MCW rank 5 bound to socket 1[core 15[hwt 0-1]]: [../../../../../../../../../../../../../..][../BB/../../../../../../../../../../../..]
[auricle:14958] MCW rank 6 bound to socket 1[core 16[hwt 0-1]]: [../../../../../../../../../../../../../..][../../BB/../../../../../../../../../../..]
[auricle:14958] MCW rank 7 bound to socket 1[core 17[hwt 0-1]]: [../../../../../../../../../../../../../..][../../../BB/../../../../../../../../../..]
[auricle:14958] MCW rank 0 bound to socket 0[core 0[hwt 0-1]]:  [BB/../../../../../../../../../../../../..][../../../../../../../../../../../../../..]
[auricle:14958] MCW rank 1 bound to socket 0[core 1[hwt 0-1]]:  [../BB/../../../../../../../../../../../..][../../../../../../../../../../../../../..]
[auricle:14958] MCW rank 2 bound to socket 0[core 2[hwt 0-1]]:  [../../BB/../../../../../../../../../../..][../../../../../../../../../../../../../..]
[auricle:14958] MCW rank 3 bound to socket 0[core 3[hwt 0-1]]:  [../../../BB/../../../../../../../../../..][../../../../../../../../../../../../../..]
8  48484.6851   Rate (MB/s) 4.05238
```
Note that this defaults to using all known sockets automatically, but if we
supply `-n 4` as well then we obtain
```
$ mpiexec -n 4 --bind-to core -report-bindings --map-by ppr:4:socket ./streams
[auricle:15285] MCW rank 2 bound to socket 0[core 2[hwt 0-1]]: [../../BB/../../../../../../../../../../..][../../../../../../../../../../../../../..]
[auricle:15285] MCW rank 3 bound to socket 0[core 3[hwt 0-1]]: [../../../BB/../../../../../../../../../..][../../../../../../../../../../../../../..]
[auricle:15285] MCW rank 0 bound to socket 0[core 0[hwt 0-1]]: [BB/../../../../../../../../../../../../..][../../../../../../../../../../../../../..]
[auricle:15285] MCW rank 1 bound to socket 0[core 1[hwt 0-1]]: [../BB/../../../../../../../../../../../..][../../../../../../../../../../../../../..]
4  24219.5586   Rate (MB/s) 2.02429
```

i.e., we get a single socket. If we ask it to do something impossible we get a
useful error message:
```
$ mpiexec -n 8 --bind-to core -report-bindings --map-by ppr:2:socket ./streams
--------------------------------------------------------------------------
Your job has requested more processes than the ppr for
this topology can support:

  App: ./streams
  Number of procs:  8
  PPR: 2:socket

Please revise the conflict and try again.
--------------------------------------------------------------------------
```

### Setting up scripts for SLURM
This section assumes that you are using SLURM to manage running IBAMR on a
cluster.

Scientific programs, like IBAMR applications, typically run for long periods of
time and require known amount of resources: for example, one might want to run
their IBAMR application with 100 processors for four days.

To allocate resources fairly most clusters use a scheduling program to assign
jobs to actual computational hardware.

Resource allocation typically require submitting a shell script describing both
what is required (e.g., amount of time, number of nodes, etc.) and also a
description of the work. Here is an example job script for running streams:
```
#!/bin/bash
#SBATCH --partition=debug_queue
#SBATCH --nodes=2
#SBATCH --time=00:10:00
#SBATCH --job-name=run-streams
# When to send an email - here, send an email when the job starts, finishes,
# or encounters an error.
#SBATCH --mail-type=begin,end,error
# Email address to use.
#SBATCH --email-user=user@email.com
# not necessary for most queues on dogwood (those default to exclusive node
# access), but it can't hurt to ask. This guarantees that we are the only
# one running things on the provided nodes. Since we typically do not use
# every processor (but we do want access to all the memory bandwidth) its
# important to provide this on machines where it is not the default.
#SBATCH --exclusive

module load openmpi_4.0.1/gcc_9.1.0

cd ~/Documents/Code/C/streams
mpiexec -map-by ppr:10:socket -bind-to core ./streams
```
Here `mpiexec` knows to work with SLURM to create a total of 40 processes -
i.e., it will automatically detect the number of nodes and start twenty
processes on each node (and ten on each socket).

You can submit your job to the job scheduler by running `sbatch script.sb`,
where `script.sb` is the file name of the script you want to run (such as the
script given above).

### SLURM tips and tricks
#### keeping track of jobs
Try running `squeue -u username` (where `username` is your username) to see the
current status of your jobs. Running
```
watch -n 1 squeue -u username
```
will update the result every second.

#### viewing all jobs
Try running `squeue` to print *all* current jobs on the cluster. Try running
`finger username` to get more information on a certain user.

### Final advice
From the PETSc manual:

"For a typical, memory bandwidth-limited PETSc application, the primary
consideration in placing MPI processes is ensuring that processes are evenly
distributed among sockets, and hence using all available memory channels."

The same holds for IBAMR: one should use the streams benchmark to determine how
many processes should be put on a socket and then use `--map-by ppr:N:socket` to
set up a run with `N` processes per socket. `--bind-to core` is a good default
but might not be the optimal choice on all hardware.
