#!/bin/bash
# mpnn/run_mpnn.sh
set -euo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: $0 <input_files_timestamp_folder>"
  exit 1
fi

INPUT_DIR="$(cd "$1" && pwd)"

source ~/anaconda3/etc/profile.d/conda.sh
conda activate ligandmpnn_env

OUT_DIR="$(dirname "$(dirname "$INPUT_DIR")")/mpnn_out/$(basename "$INPUT_DIR")"
mkdir -p "$OUT_DIR"

cd /home/panda/LigandMPNN #TODO: Point to LigandMPNN repository.

python run.py \
  --temperature 0.15 \
  --verbose 0 \
  --pdb_path_multi "$INPUT_DIR/pdb_ids.json" \
  --redesigned_residues_multi "$INPUT_DIR/redesigned_residues_multi.json" \
  --out_folder "$OUT_DIR" \
  --seed 42 \
  --batch_size 10 \
  --number_of_batches 5 \
  --file_ending "_a_et_b"
