#!/usr/bin/env python3
"""Generate AF3 v4 homomer inputs for FHV protein A."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


CHAIN_IDS = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")


def read_fasta(path: Path) -> str:
    lines = [line.strip() for line in path.read_text().splitlines()]
    sequence = "".join(line for line in lines if line and not line.startswith(">"))
    invalid = sorted(set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"))
    if invalid:
        raise ValueError(f"Unsupported FASTA characters: {invalid}")
    return sequence


def write_input(out_dir: Path, sequence: str, copies: int, seeds: list[int]) -> Path:
    if not 1 <= copies <= len(CHAIN_IDS):
        raise ValueError(f"copies must be between 1 and {len(CHAIN_IDS)}")
    label = {1: "monomer", 2: "dimer", 4: "tetramer"}.get(copies, f"c{copies}")
    name = f"fhv_a_{label}"
    entity_id: str | list[str]
    entity_id = CHAIN_IDS[0] if copies == 1 else CHAIN_IDS[:copies]
    payload = {
        "name": name,
        "modelSeeds": seeds,
        "sequences": [
            {
                "protein": {
                    "id": entity_id,
                    "sequence": sequence,
                    "description": f"FHV protein A full-length homomer; {copies} copies",
                }
            }
        ],
        "dialect": "alphafold3",
        "version": 4,
    }
    path = out_dir / f"{name}.json"
    path.write_text(json.dumps(payload, indent=2) + "\n")
    return path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, default=Path("reference/fhv_protein_a.fasta"))
    parser.add_argument("--out-dir", type=Path, default=Path("inputs"))
    parser.add_argument("--copies", type=int, nargs="+", default=[1, 2, 4])
    parser.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    args = parser.parse_args()

    sequence = read_fasta(args.fasta)
    if len(sequence) != 998:
        raise ValueError(f"Expected the 998-aa 8FM9 protein A sequence, got {len(sequence)} aa")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for copies in args.copies:
        tokens = copies * len(sequence)
        if tokens > 5120:
            print(
                f"WARNING: {copies} copies = {tokens} protein tokens, above AF3's "
                "published 5,120-token A100/H100 80 GB validation size."
            )
        print(write_input(args.out_dir, sequence, copies, args.seeds))


if __name__ == "__main__":
    main()

