#!/usr/bin/env bash
set -euo pipefail

job_name=${1:?usage: run_af3_inference.sh JOB_NAME [INPUT_JSON]}
input_override=${2:-}
: "${AF3_ROOT:?AF3_ROOT must be set by the submit file}"

model_dir="$AF3_ROOT/models"
output_root=${AF3_OUTPUT_ROOT:-"$AF3_ROOT/runs"}
job_output_dir="$output_root/$job_name"
input_json="$job_output_dir/${job_name}_data.json"
cache_dir="$AF3_ROOT/jax-cache"

if [[ -n "$input_override" ]]; then
  input_json="$input_override"
fi

# Compatibility with earlier wrappers that passed a job-specific output
# directory to AF3, causing AF3 to create a second nested job directory.
input_candidates=(
  "$job_output_dir/$job_name/${job_name}_data.json"
  "$AF3_ROOT/runs/$job_name/${job_name}_data.json"
  "$AF3_ROOT/runs/$job_name/$job_name/${job_name}_data.json"
)
if [[ ! -s "$input_json" ]]; then
  for candidate in "${input_candidates[@]}"; do
    if [[ -s "$candidate" ]]; then
      input_json="$candidate"
      break
    fi
  done
fi
[[ -s "$input_json" ]] || {
  echo "ERROR: processed AF3 input is missing or empty: $input_json" >&2
  exit 2
}
[[ -s "$model_dir/af3.bin.zst" ]] || {
  echo "ERROR: AF3 model weights are missing: $model_dir/af3.bin.zst" >&2
  echo "The database is not needed for this inference-only job." >&2
  exit 2
}
mkdir -p "$output_root" "$cache_dir"
echo "AF3 processed input: $input_json"
echo "AF3 prediction output: $job_output_dir"
nvidia-smi

python /app/alphafold/run_alphafold.py \
  --json_path="$input_json" \
  --model_dir="$model_dir" \
  --output_dir="$output_root" \
  --norun_data_pipeline \
  --force_output_dir \
  --jax_compilation_cache_dir="$cache_dir" \
  --compress_large_output_files

if [[ -s "$job_output_dir/${job_name}_model.cif.zst" ]]; then
  model_output="$job_output_dir/${job_name}_model.cif.zst"
elif [[ -s "$job_output_dir/${job_name}_model.cif" ]]; then
  model_output="$job_output_dir/${job_name}_model.cif"
else
  echo "ERROR: AF3 finished without a top-ranked model CIF in $job_output_dir" >&2
  exit 2
fi
echo "AF3 top-ranked model: $model_output"
