#!/bin/bash
#SBATCH -p biodroid          # Partition name
#SBATCH -c 10                 # Number of CPU cores
#SBATCH --gres=gpu:1         # Request one GPU
#SBATCH --mem=46G            # Memory allocation
#SBATCH --time=10-00:00:00      # Time limit

module load biodroid
module load tmalign/1.0
module load mmseqs2/1.0
eval "$(conda shell.bash hook)"
conda activate /projects/biodroid/conda_envs/miniconda3/envs/mpnn_retrain

python3 ../training.py \
    --path_for_training_data ../pdb_pt \
    --path_for_data_files ../data_files \
    --path_for_outputs ../010_train_output \
    --num_neighbors 48 \
    --backbone_noise 0.10 \
    --num_epochs 300 \
    --batch_size 30000 \
    --max_protein_length 8000 \
    --num_examples_per_epoch 10000 \
    --hidden_dim 128 \
    --suffix cryo_010_model \
    --rerun


