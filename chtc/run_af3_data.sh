#!/usr/bin/env bash
set -euo pipefail

input_json=${1:?usage: run_af3_data.sh INPUT_JSON JOB_NAME}
job_name=${2:?usage: run_af3_data.sh INPUT_JSON JOB_NAME}
: "${AF3_ROOT:?AF3_ROOT must be set by the submit file}"

db_dir="$AF3_ROOT/databases"
output_dir="$AF3_ROOT/runs/$job_name"
mkdir -p "$output_dir"

python /app/alphafold/run_alphafold.py \
  --json_path="$input_json" \
  --db_dir="$db_dir" \
  --output_dir="$output_dir" \
  --norun_inference \
  --force_output_dir

test -s "$output_dir/${job_name}_data.json"

