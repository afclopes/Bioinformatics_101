# How is the cluster organised?

There are 3 types of folders:
 1. Raw folders: to save your fastq data that is not processed
 2. Work folders: to save temporary files that you are working on
 3. Data folders: to save final processed output files which must be backed up
 
 To check the amount of space leftover in each directory: go to the directory, and do >> df -h .

# Which type of partitions are there on the cluster?

 - CPU partitions only
 
 - CPU and GPU partitions: GPU-accelerated non-interactive computations
 
 - GPU partition: includes the "CPU and GPU partition" nodes, but also includes nodes that are actually workstations (eg. CPU nodes). 
 This partition is used for interactive, graphic-heavy computations, *OR* GPU-accelerated calculations that *have* to be launched from 
 a GUI that itself needs GPU-acceleration. 
 
 If you don't need to sit in front of a GPU-enabled workstation, then I'd avoid using the "gpu" partition completely. 
 (If you do need them, then there are ways to request a specific node rather than using a partition name).
 
 ### What are GPUs? 
 
GPU: it's a special type of CPU than is good at running many parallel (simultaneous) operations. It demands special codes. 
For example, awk code cannot run on GPUs.

A node with a GPU also has a CPU. If you submit your pipeline full of awk and gzip to a GPU node, the code will run on the CPU, but the GPU will be wasted.

People can write python scripts that uses special libraries that are run on a GPU. Then, they can have python running on the CPU but using 
libraries that run on the GPU. Normally though, if you dont use special libraries that are run on GPU, then you can use python on the CPU without any issues.

Summary: for now, stick to using the CPUs. 
