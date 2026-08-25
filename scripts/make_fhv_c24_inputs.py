#!/usr/bin/env python3
"""Generate staged AlphaFold 3 inputs for the FHV mature C24 crown.

AlphaFold 3 cannot consume an EM map.  These jobs sample tractable sequence
and local-interface hypotheses; their output must subsequently be fitted to
and refined against the experimental density.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path


AA = set("ACDEFGHIKLMNPQRSTVWY")
CHAIN_IDS = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
PUBLISHED_80GB_LIMIT = 5120


@dataclass(frozen=True)
class Experiment:
    name: str
    start: int
    end: int
    copies: int
    tier: str
    purpose: str

    @property
    def residues_per_chain(self) -> int:
        return self.end - self.start + 1

    @property
    def protein_tokens(self) -> int:
        return self.residues_per_chain * self.copies


RECOMMENDED = (
    Experiment("fhv_g4wwb0_full_monomer", 1, 998, 1, "recommended", "single-chain baseline"),
    Experiment("fhv_g4wwb0_full_dimer", 1, 998, 2, "recommended", "adjacent-interface hypotheses"),
    Experiment("fhv_g4wwb0_full_tetramer", 1, 998, 4, "recommended", "local crown neighborhood"),
    Experiment("fhv_g4wwb0_floor_core_55_396_c12", 55, 396, 12, "recommended", "C12 floor-core hypothesis"),
    Experiment("fhv_g4wwb0_nterm_1_452_c4", 1, 452, 4, "recommended", "local floor/leg alternatives"),
    Experiment("fhv_g4wwb0_pol_453_896_c4", 453, 896, 4, "recommended", "local polymerase-ring alternatives"),
    Experiment("fhv_g4wwb0_linker_pol_379_896_c4", 379, 896, 4, "recommended", "linker-to-polymerase context"),
)

RISKY_C12 = (
    Experiment("fhv_g4wwb0_nterm_1_452_c12", 1, 452, 12, "oversize", "complete C12 N-terminal ring"),
    Experiment("fhv_g4wwb0_pol_453_896_c12", 453, 896, 12, "oversize", "complete C12 polymerase ring"),
)

IMPRACTICAL_C24 = (
    Experiment("fhv_g4wwb0_full_c24", 1, 998, 24, "impractical", "direct mature-crown hypothesis"),
)


def read_fasta(path: Path) -> tuple[str, str]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"{path} is not a FASTA file")
    sequence = "".join(line.upper() for line in lines[1:] if not line.startswith(">"))
    invalid = sorted(set(sequence) - AA)
    if invalid:
        raise ValueError(f"Unsupported FASTA characters: {invalid}")
    return lines[0][1:], sequence


def make_payload(experiment: Experiment, sequence: str, seeds: list[int]) -> dict:
    fragment = sequence[experiment.start - 1 : experiment.end]
    ids: str | list[str]
    ids = CHAIN_IDS[0] if experiment.copies == 1 else CHAIN_IDS[: experiment.copies]
    copy_word = "copy" if experiment.copies == 1 else "copies"
    return {
        "name": experiment.name,
        "modelSeeds": seeds,
        "sequences": [
            {
                "protein": {
                    "id": ids,
                    "sequence": fragment,
                    "description": (
                        f"FHV protein A G4WWB0 residues {experiment.start}-{experiment.end}; "
                        f"{experiment.copies} identical {copy_word}; {experiment.purpose}"
                    ),
                }
            }
        ],
        "dialect": "alphafold3",
        "version": 4,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, default=Path("reference/FHV_proteinA_G4WWB0.fasta"))
    parser.add_argument("--out-dir", type=Path, default=Path("inputs/fhv_g4wwb0"))
    parser.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    parser.add_argument(
        "--include-risky-c12",
        action="store_true",
        help="also write the 5,328- and 5,424-token C12 domain jobs",
    )
    parser.add_argument(
        "--include-impractical-c24",
        action="store_true",
        help="write the 23,952-token full C24 JSON (not recommended for CHTC)",
    )
    args = parser.parse_args()

    header, sequence = read_fasta(args.fasta)
    if len(sequence) != 998:
        raise ValueError(f"Expected 998-aa FHV protein A, got {len(sequence)} aa")
    if not args.seeds or any(seed < 0 for seed in args.seeds):
        raise ValueError("Provide one or more non-negative integer seeds")

    experiments = list(RECOMMENDED)
    if args.include_risky_c12:
        experiments.extend(RISKY_C12)
    if args.include_impractical_c24:
        experiments.extend(IMPRACTICAL_C24)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for experiment in experiments:
        if experiment.end > len(sequence):
            raise ValueError(f"{experiment.name} ends after residue {len(sequence)}")
        payload = make_payload(experiment, sequence, args.seeds)
        output = args.out_dir / f"{experiment.name}.json"
        output.write_text(json.dumps(payload, indent=2) + "\n")
        status = "within published 80-GB benchmark" if experiment.protein_tokens <= PUBLISHED_80GB_LIMIT else "oversize"
        rows.append(
            {
                "name": experiment.name,
                "residue_range": f"{experiment.start}-{experiment.end}",
                "copies": experiment.copies,
                "protein_tokens": experiment.protein_tokens,
                "tier": experiment.tier,
                "status": status,
                "purpose": experiment.purpose,
                "json": output.name,
            }
        )
        print(f"{output}: {experiment.protein_tokens} protein tokens ({status})")

    manifest = args.out_dir / "manifest.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    metadata = args.out_dir / "source.json"
    metadata.write_text(
        json.dumps(
            {
                "fasta": str(args.fasta),
                "header": header,
                "sequence_length": len(sequence),
                "seeds": args.seeds,
                "note": "Residue ranges are inclusive and refer to the supplied 998-aa G4WWB0 sequence.",
            },
            indent=2,
        )
        + "\n"
    )
    print(manifest)


if __name__ == "__main__":
    main()
