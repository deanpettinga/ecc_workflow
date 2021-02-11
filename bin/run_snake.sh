#!/bin/bash
#SBATCH -J circle
#SBATCH --partition=savio
#SBATCH --account=fc_kvkallow
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH -o logs/runs/circle.%j.o
#SBATCH -e logs/runs/circle.%j.e
#SBATCH --mail-type=All

source ~/.bashrc
conda activate snakemake

# save DAG job file with time stamp
TIME=$(date "+%Y-%m-%d_%H.%M.%S")

# make logs dir if it does not exist already. Without this, logs/ is automatically generate only after the first run of the pipeline
logs_dir="logs/runs"
[[ -d $logs_dir ]] || mkdir -p $logs_dir

snakemake --use-conda -n > logs/runs/workflow_${TIME}.txt
snakemake --dag | dot -Tpng > logs/runs/workflow_${TIME}.png

snakemake \
--use-conda \
--jobs 20 \
--cluster "sbatch \
  --partition=savio \
  --account=fc_kvkallow \
  --qos=savio_normal \
  --nodes={resources.nodes} \
  --cpus-per-task={resources.threads} \
  --time=72:00:00"
