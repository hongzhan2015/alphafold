#!/usr/bin/env bash
set -euo pipefail

input_json=${1:?usage: run_af3_data.sh INPUT_JSON JOB_NAME}
job_name=${2:?usage: run_af3_data.sh INPUT_JSON JOB_NAME}
: "${AF3_ROOT:?AF3_ROOT must be set by the submit file}"

# Prefer CHTC's full AF3 database installation on nodes advertising
# HasAlphafold3. Fall back to a user-managed database only when present.
if [[ -n "${AF3_DB_DIR:-}" ]]; then
  db_dir="$AF3_DB_DIR"
elif [[ -f /alphafold3/uniref90_2022_05.fa ]]; then
  db_dir=/alphafold3
elif [[ -f "$AF3_ROOT/databases/uniref90_2022_05.fa" ]]; then
  db_dir="$AF3_ROOT/databases"
elif [[ -f "$AF3_ROOT/public_databases/uniref90_2022_05.fa" ]]; then
  db_dir="$AF3_ROOT/public_databases"
else
  echo "ERROR: no AlphaFold 3 database installation was found." >&2
  echo "Expected CHTC shared databases under /alphafold3 on a HasAlphafold3 node." >&2
  exit 2
fi
output_root=${AF3_OUTPUT_ROOT:-"$AF3_ROOT/runs"}
output_dir="$output_root/$job_name"
[[ -s "$input_json" ]] || {
  echo "ERROR: transferred AF3 input is missing or empty: $input_json" >&2
  echo "Working directory: $(pwd)" >&2
  exit 2
}
mkdir -p "$output_dir"
echo "AF3 database directory: $db_dir"

# The sequence is carried in the AF3 JSON rather than embedded in this script.
# Print a concise check in the Condor output before starting the long pipeline.
python - "$input_json" "$job_name" <<'PY'
import json
import sys

path, expected_name = sys.argv[1:]
with open(path) as handle:
    payload = json.load(handle)
if payload.get("name") != expected_name:
    raise SystemExit(
        f"JSON name {payload.get('name')!r} does not match job_name {expected_name!r}"
    )
proteins = [entry["protein"] for entry in payload.get("sequences", []) if "protein" in entry]
if not proteins:
    raise SystemExit("The input JSON has no protein sequence")
for protein in proteins:
    chain_ids = protein["id"] if isinstance(protein["id"], list) else [protein["id"]]
    print(
        f"AF3 input: {expected_name}; chains={','.join(chain_ids)}; "
        f"residues_per_chain={len(protein['sequence'])}; "
        f"protein_tokens={len(protein['sequence']) * len(chain_ids)}",
        flush=True,
    )
PY

python /app/alphafold/run_alphafold.py \
  --json_path="$input_json" \
  --db_dir="$db_dir" \
  --output_dir="$output_dir" \
  --norun_inference \
  --force_output_dir

test -s "$output_dir/${job_name}_data.json"
echo "AF3 processed input: $output_dir/${job_name}_data.json"
