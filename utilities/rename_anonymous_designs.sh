#!/usr/bin/env bash
cd /projects/arvid/c_symmetry_design/rfd2/out || exit 1

for dir in */; do
  config_id="${dir%/}"
  for f in "$dir"design_*.pdb; do
    [ -e "$f" ] || continue
    mv "$f" "${dir}${config_id}_$(basename "$f")"
  done
  for f in "$dir"design_*.trb; do
    [ -e "$f" ] || continue
    mv "$f" "${dir}${config_id}_$(basename "$f")"
  done
done
