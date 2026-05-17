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

mkdir -p ../cluster_data

# create fasta file
python3 -W ignore ../scripts/get_multifasta_from_pdb_path.py \
	-p ../input_pdb/ -c A -o ../cluster_data/input_pdb_seqs.fasta

# Clustering the fasta sequences
OMP_PROC_BIND=false
echo "Start clustering"
mmseqs easy-cluster ../cluster_data/input_pdb_seqs.fasta ../cluster_data/input_pdb_seqs ./tmp \
 	--min-seq-id 0.3 -c 0.8 --cov-mode 0
