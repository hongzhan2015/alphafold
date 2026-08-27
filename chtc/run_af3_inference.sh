#!/usr/bin/env bash
set -euo pipefail

job_name=${1:?usage: run_af3_inference.sh JOB_NAME}
: "${AF3_ROOT:?AF3_ROOT must be set by the submit file}"

model_dir="$AF3_ROOT/models"
output_root=${AF3_OUTPUT_ROOT:-"$AF3_ROOT/runs"}
output_dir="$output_root/$job_name"
input_json="$output_dir/${job_name}_data.json"
cache_dir="$AF3_ROOT/jax-cache"

# Compatibility with data-pipeline jobs submitted before AF3_OUTPUT_ROOT was
# introduced: read their processed JSON from the previous runs directory while
# writing all inference products to the requested output directory.
if [[ ! -s "$input_json" && -s "$AF3_ROOT/runs/$job_name/${job_name}_data.json" ]]; then
  input_json="$AF3_ROOT/runs/$job_name/${job_name}_data.json"
fi
test -s "$input_json"
mkdir -p "$output_dir" "$cache_dir"
echo "AF3 processed input: $input_json"
echo "AF3 prediction output: $output_dir"
nvidia-smi

python /app/alphafold/run_alphafold.py \
  --json_path="$input_json" \
  --model_dir="$model_dir" \
  --output_dir="$output_dir" \
  --norun_data_pipeline \
  --force_output_dir \
  --jax_compilation_cache_dir="$cache_dir" \
  --compress_large_output_files
