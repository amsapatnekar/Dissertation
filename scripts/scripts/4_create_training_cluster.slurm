#!/bin/bash
#SBATCH -p biodroid          # Partition name
#SBATCH -c 10                 # Number of CPU cores
#SBATCH --gres=gpu:1         # Request one GPU
#SBATCH --mem=46G            # Memory allocation
#SBATCH --time=03-00:00:00      # Time limit

module load biodroid
module load tmalign/1.0
module load mmseqs2/1.0
eval "$(conda shell.bash hook)"
conda activate /projects/biodroid/conda_envs/miniconda3/envs/mpnn_retrain

python3 ./create_cluster_lists.py ../pdb_pt/ ../cluster_data/input_pdb_seqs_cluster.tsv
