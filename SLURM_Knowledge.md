# Must knowns of SLURM: running jobs on the cluster

Extra knowledge can be found here: https://slurm.schedmd.com/sbatch.html

1. Essential commands
    - sinfo: to tell which nodes are there
    - squeue: to tell the state of jobs running
    - sbatch: to prepare a job to run
    - salloc: to allocate CPUs
    - scancel numberjob: cancels the job with the number specified
    
What is a job: any script, code in bash, C, python that you need to run.
    
With *sbatch* you can create a job by stating which script to run and how many CPUs to use. You can then keep checking when your job is done by using *squeue*.
With *salloc* you will be asking for the available CPUs and then you will use *srun* in a live interactive session.

2. *sinfo* gives you information about which partitions or nodes are available. What to remember?
    - AVAIL: UP: means available, DOWN: means not available
    - STATE: *down\**: means not responding, *mix*: some nodes free and some nodes not free. *alloc*: allocated for a job. *idle*: free to use. *drain*: under 
    construction/fixing.
    - NODELIST: adev[1-2]: means whatever the state is, this will refer to the nodes in square brackets. so STATE down, NODELIST adev[1-2] means that the nodes 1 
    and 2 are down.
    - PARTITION: clustername\*: means that some of the nodes in that cluster are in the state of *drain* but you can still use this cluster and will be directed 
    to other nodes. Just use the clustername without the \*
    
3. *squeue* gives information about which jobs are being run, and queued. What to remember?
    - ST: state of the job. R: means running, PD: means pending
    - TIME: how long the job has already run for
    - NODELIST(REASON): which nodes are being used for the job, or what is the reason for the job remaining pending. Usual reasons are Resources 
    (waiting for resources to become available) and Priority (queued behind a higher priority job)
    - -u username : to check if there are any jobs being run by a particular user
    - -t status: eg. -t RUNNING: to see which jobs are running
   
4. *scontrol* to find information about the job you submitted in detail. What to remember?
    - scontrol show jobid -dd jobidnumber: gives you detailed information of the job and the exact script that you run (only if the job is still running) so 
    you can see what jobs are left on the cluster running still.

5. Nodes versus core (or cpu):
    
    Every node is made up of 16, 32, 56, or 64 cores. You can use a maximum of 4 nodes on this cluster. Which means you can use a very high number of cores. 
    Usually, by just adding "#SBATCH -c 4" to your script, you will already have your job running very fast. To check the number of nodes and cores available per 
    partition: > sinfo -Nl
    
6. Running a job with job-array:
    
    You can use --array to submit a job that can run in parallel multiple samples at the same time. For example:
    > sbatch --array1-12%2 filenameScript.sh
    
    Where '1-10' refers to the length of the file of sample names. And '%2' refers to how many jobs should be run in parallel. Sbatch array can be useful 
    in that you will find out which sample had an issue. But, it can be much faster to run multiple samples without the sbatch array (muuuuch faster; using 
    xargs it took a weekend to run 12 samples, starting on Friday, but when using --array it took 1 week from Friday to Friday to complete), for example by 
    using inside your sbatch script:
    
    >cat Multiple_Oocyte_List.dat | xargs -i -t -n 1 -P 4 /path/SBatch-full-align-pipeline_Multiple_Oocytes.sh {}
    
    Here, 'xargs' takes standard input and executes it for example in another script. 
    
    -i {}: Places occurrences of sample names read by 'cat' line by line from standard input. 
    
    -t: Print the command line on the standard error output before executing it. 
    
    -P: number of processes that can be run at a time (here 4 samples were run through the script at the same time). -P needs to be used always with -n.
    
    For more details on xargs: https://man7.org/linux/man-pages/man1/xargs.1.html
    
 7. Checking old jobs run on cluster:
 
 > sacct --user=afcl2 --starttime=2020-12-04 --format=jobid,jobname%50,elapsed,ntasks,exitcode,state

'sacct' allows you to look at all the jobs that are associated to that username but only for the current day. To look on previous dates, you need to specify 
the start date for the job. A selected number of columns/information will appear by default, but these often are not the most informative insights. For more 
insights, you can use '--long', but often too many columns of information usually appears. In that case, its best to use '--format' to select exactly which piece 
of information you are interested in. For more information of the options: https://curc.readthedocs.io/en/latest/running-jobs/slurm-commands.html

'jobname\%50' refers to the number of characters that you can see for the job name. This can be super useful if you job name is very long and you cannot 
distinguish the truncated names that often appear with 'squeue'.
