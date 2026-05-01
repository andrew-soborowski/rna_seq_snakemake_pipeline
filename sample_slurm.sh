#!/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=andrew.soborowski@duke.edu
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH --account=schmidlab
#SBATCH --partition=schmidlab
#SBATCH --time=48:00:00


#Be sure to change your source and conda to your own conda shell script and your snakemake environment.
#Note that this is designed to run on account=schmidlab, so if you want to run on general you also have to edit profile/config.yaml to change the account there as well.
source /hpc/group/schmidlab/als185/miniforge3/etc/profile.d/mamba.sh
mamba activate /hpc/group/schmidlab/als185/miniforge3/envs/snakemake_env/

#Recommend to do a dry run first, add a -n or --dry-run flag 
snakemake --use-conda --profile ./profile -p -k -w  30 --rerun-incomplete &> snakemake_run.log
