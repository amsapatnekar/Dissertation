#!/bin/bash
#SBATCH -p biodroid          # Partition name
#SBATCH -c 10                 # Number of CPU cores
#SBATCH --gres=tmp:32G,gpu:1         # Request one GPU
#SBATCH --mem=46G            # Memory allocation
#SBATCH --time=10-00:00:00      # Time limit
#SBATCH --array=0-242%4
#SBATCH --nice=500 
module load biodroid
module load tmalign/1.0
module load mmseqs2/1.0
eval "$(conda shell.bash hook)"
conda activate /projects/biodroid/conda_envs/miniconda3/envs/mpnn_retrain

num_examples_per_epoch=(15000 20000 30000)
batch_size=(8500 10000 15000)
num_neighbours=(32 48 58)
hidden_dim=(64 128 156)
dropout=(0.05 0.1 0.15)

i=0
for bs in "${batch_size[@]}"; do
    for ex in "${num_examples_per_epoch[@]}"; do
        for nn in "${num_neighbours[@]}"; do
            for hd in "${hidden_dim[@]}"; do
                for dp in "${dropout[@]}"; do
                    if [ "$i" = "$SLURM_ARRAY_TASK_ID" ]; then
                        BATCH_SIZE=$bs
                        EXAMPLES=$ex
                        NEIGHBOURS=$nn
                        HIDDEN_DIM=$hd
                        DROPOUT=$dp
                    fi
                    ((i++))
                done
            done
        done
    done
done

SUFFIX="gridsearch_${SLURM_ARRAY_TASK_ID}"
echo "Job $SLURM_ARRAY_TASK_ID running with:"
echo "  Batch Size:      $BATCH_SIZE"
echo "  Examples/Epoch:  $EXAMPLES"
echo "  Neighbors:       $NEIGHBOURS"
echo "  Hidden Dim:      $HIDDEN_DIM"
echo "  Dropout:         $DROPOUT"
echo "  Suffix:          $SUFFIX"

python3 ../training.py \
    --path_for_training_data ../pdb_pt \
    --path_for_data_files ../data_files \
    --path_for_outputs ../train_output \
    --num_neighbors $NEIGHBOURS \
    --backbone_noise 0.20 \
    --num_epochs 300 \
    --batch_size $BATCH_SIZE \
    --max_protein_length 8000 \
    --num_examples_per_epoch $EXAMPLES \
    --suffix $SUFFIX \
    --hidden_dim $HIDDEN_DIM \
    --dropout $DROPOUT \
    --rerun
