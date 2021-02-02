# Anaconda

When using anaconda, it is worth keeping in mind that there are multiple variantions of this, eg. anaconda3, miniconda2. When installing or updating things, you must ensure that you are in the correct path, since the commands for this are the same for all different versions of 'conda'.

It is worth noting that sometimes it can be quite useful to create a new environment where you can install specific packages that you believe you will need.

Example: 

> conda create xxx pythona.b --> 'pythona.b' means the version of python you want. For example 3.8 or 3.9 or 3.6 'xxx' is the name of you new environment

> conda activate py3.8 --> activate new environment

> echo $PATH --> check that this new environment is in your path

> conda install -c bioconda bcftools --> install your specific package

> which bcftools --> check that the package got installed in your path

> bcftools --version --> check the version of your package is up to date

> conda list -n myenv --> check which packages are installed in my environment 'myenv'

Interestingly, the moment that you do 'conda activate' on your area, and submit an sbatch job, then the server will copy everything that is in your path before it runs your jobs.


To find out which packages you can download from anaconda: https://anaconda.org/search?q=bowtie2
