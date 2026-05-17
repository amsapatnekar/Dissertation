#!/bin/bash
#SBATCH -p biodroid          # Partition name
#SBATCH -c 10                 # Number of CPU cores
#SBATCH --gres=gpu:1         # Request one GPU
#SBATCH --mem=46G            # Memory allocation
#SBATCH --time=03-00:00:00      # Time limit

module load biodroid
module load tmalign/1.0
eval "$(conda shell.bash hook)"
conda activate /projects/biodroid/conda_envs/miniconda3/envs/mpnn_retrain

n_splits=2000
n_process=20

mkdir -p ../chunk
split -l $n_splits ../data_files/id_map.txt ../chunk/chunk_
# ONLY RUN ONES !!!
ls ../chunk/chunk_* | xargs -n 1 -P $n_process bash -c 'while read file; do bash helper_create_pt_files.sh $file; done < $1' _
