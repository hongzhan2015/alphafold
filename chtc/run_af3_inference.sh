#!/usr/bin/env bash
set -euo pipefail

job_name=${1:?usage: run_af3_inference.sh JOB_NAME}
: "${AF3_ROOT:?AF3_ROOT must be set by the submit file}"

model_dir="$AF3_ROOT/models"
output_dir="$AF3_ROOT/runs/$job_name"
input_json="$output_dir/${job_name}_data.json"
cache_dir="$AF3_ROOT/jax-cache"

test -s "$input_json"
mkdir -p "$output_dir" "$cache_dir"
nvidia-smi

python /app/alphafold/run_alphafold.py \
  --json_path="$input_json" \
  --model_dir="$model_dir" \
  --output_dir="$output_dir" \
  --norun_data_pipeline \
  --force_output_dir \
  --jax_compilation_cache_dir="$cache_dir" \
  --compress_large_output_files

