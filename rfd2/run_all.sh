#!/bin/bash
# rfd2/run_all.sh
set -o pipefail

set +u
source ~/anaconda3/etc/profile.d/conda.sh
conda activate ca_rfd
set -u

export PYTHONPATH=/home/panda/Resources/software/rfd2

CONFIG_DIR="$(pwd)/configs/generated"
LOG_DIR="$(pwd)/logs"
RFD2_SCRIPT="/home/panda/Resources/software/rfd2/rf_diffusion/run_inference.py"
mkdir -p "$LOG_DIR"

cd /home/panda/Resources/software/rfd2

for config in "$CONFIG_DIR"/*.yaml; do
    name=$(basename "$config" .yaml)
    echo "Running $name..."
    if python "$RFD2_SCRIPT" --config-dir "$CONFIG_DIR" --config-name "$name" \
        > "$LOG_DIR/${name}.log" 2>&1; then
        echo "  OK: $name"
    else
        echo "  FAILED: $name (see $LOG_DIR/${name}.log)"
    fi
done
