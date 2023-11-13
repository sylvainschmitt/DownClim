#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH -J smkTemplate
#SBATCH -o smkTemplate.%N.%j.out
#SBATCH -e smkTemplate.%N.%j.err
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --account=agap
#SBATCH --partition=agap_normal #sinfo

# Environment
module purge
module load snakemake
module load singularity

# Variables
CONFIG=config/ressources.muse.yaml
COMMAND="sbatch --cpus-per-task={cluster.cpus} --time={cluster.time} --mem={cluster.mem} -J {cluster.jobname} -o snake_subjob_log/{cluster.jobname}.%N.%j.out -e snake_subjob_log/{cluster.jobname}.%N.%j.err --account=agap --partition={cluster.partition}"
CORES=100
mkdir -p snake_subjob_log

# Workflow
snakemake -s Snakefile --use-singularity -j $CORES --cluster-config $CONFIG --cluster "$COMMAND" --keep-going # maybe --singularity-args "\-\-containall" will be need to avoid conflict with local R env

## Session informations
echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job ID:' $SLURM_JOB_ID
echo 'Number of nodes assigned to job:' $SLURM_JOB_NUM_NODES
echo 'Nodes assigned to job:' $SLURM_JOB_NODELIST
echo 'Directory:' $(pwd)
echo '########################################'
