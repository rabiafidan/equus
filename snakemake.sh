#!/bin/bash -l
#SBATCH --partition=macaque1                # write partition name(chimp;bonobo;macaque$i)
#SBATCH --job-name=snakemake                 # write your job name
#SBATCH --ntasks=1                   # number of task
#SBATCH --output=1.txt
#SBATCH --error=1.err


echo "SLURM_NODELIST $SLURM_NODELIST"
echo "NUMBER OF CORES $SLURM_NTASKS"

source ~/.bashrc

snakemake --latency-wait 300 --rerun-incomplete --cluster "sbatch --mem {resources.mem_mb} --ntasks {threads} -p macaque4 --job-name snakemake --error {params.err} --output {params.out}" --jobs 20 --keep-going



exit
