#!/usr/bin/env bash
#SBATCH -J Boltz2
#SBATCH -p biodroid
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 1-00:00:00
#SBATCH -o slurm_output/%x_%j.out
#SBATCH --export=NONE

set -euo pipefail

# ==============================
# User-configurable variables
# ==============================

BASE_DIR=$PWD
INPUT_DIR="$BASE_DIR/input_pdb"
OUT_DIR="$BASE_DIR/boltz_output"
BOLTZ_CACHE="$BASE_DIR/boltz_cache"

mkdir -p "$OUT_DIR"
mkdir -p "$BOLTZ_CACHE"
mkdir -p slurm_output

# ==============================
# Locate YAML file(s)
# ==============================

YAML_FILES=("$INPUT_DIR"/*.yaml)

if [ ${#YAML_FILES[@]} -eq 0 ] || [ ! -f "${YAML_FILES[0]}" ]; then
    echo "ERROR: No YAML files found in $INPUT_DIR"
    exit 1
fi

echo "Found ${#YAML_FILES[@]} YAML file(s) in $INPUT_DIR"

# ==============================
# Load modules
# ==============================

module purge
module load biodroid
module load boltz

# ==============================
# Run Boltz-2 on each YAML
# ==============================

for YAML_FILE in "${YAML_FILES[@]}"; do
    BASENAME=$(basename "$YAML_FILE" .yaml)
    OUT_PATH="${OUT_DIR}/${BASENAME}"
    mkdir -p "$OUT_PATH"

    echo "========================================"
    echo "[boltz] Running: $BASENAME"
    echo "[boltz] Input:   $YAML_FILE"
    echo "[boltz] Output:  $OUT_PATH"
    echo "========================================"

    boltz predict "$YAML_FILE" \
        --out_dir "$OUT_PATH" \
        --recycling_steps 15 \
        --diffusion_samples 10 \
        --cache "$BOLTZ_CACHE" \
        --num_workers 1 \
        --output_format pdb \
        --write_full_pae \
        --write_full_pde \
        --use_potentials

    echo "[boltz] Finished: $BASENAME"
done

echo "DONE: Boltz-2 completed for ${#YAML_FILES[@]} input(s)"
