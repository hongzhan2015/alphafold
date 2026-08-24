#!/usr/bin/env python3
"""Coarse rigid-body PDB-to-MRC fit for an already axis-aligned C12 assembly.

This is an initial docking utility, not a substitute for map-refinement software.
It searches Z-axis rotation, translation, and optional reversal of the model Z
direction using CA-to-density sampling.
"""

from __future__ import annotations

import argparse
import json
import math
import struct
from pathlib import Path

import numpy as np


def mrc_header(path: Path) -> tuple[tuple[int, int, int], np.ndarray, np.ndarray, int]:
    raw = path.read_bytes()[:1024]
    ints = struct.unpack("<256i", raw)
    floats = struct.unpack("<256f", raw)
    nx, ny, nz, mode = ints[:4]
    if mode != 2:
        raise ValueError("This utility currently requires MRC mode 2 (float32)")
    mx, my, mz = ints[7:10]
    cell = np.array(floats[10:13], dtype=float)
    voxel = cell / np.array([mx, my, mz], dtype=float)
    origin = np.array(floats[49:52], dtype=float)
    return (nx, ny, nz), voxel, origin, ints[23]


def pdb_ca(path: Path) -> tuple[np.ndarray, np.ndarray]:
    coords, groups = [], []
    for line in path.read_text().splitlines():
        if line.startswith("ATOM  ") and line[12:16].strip() == "CA":
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            groups.append(0 if line[21] in "ABCDEFGHIJKL" else 1)
    if not coords:
        raise ValueError("No CA atoms found")
    return np.asarray(coords, dtype=float), np.asarray(groups, dtype=int)


def transform(
    coords: np.ndarray,
    center: np.ndarray,
    target: np.ndarray,
    angle_deg: float,
    reverse_z: bool,
) -> np.ndarray:
    q = coords - center
    if reverse_z:
        q[:, 2] *= -1
    angle = math.radians(angle_deg)
    c, s = math.cos(angle), math.sin(angle)
    x = c * q[:, 0] - s * q[:, 1]
    y = s * q[:, 0] + c * q[:, 1]
    return np.column_stack((x, y, q[:, 2])) + target


def score_points(
    points: np.ndarray,
    data: np.ndarray,
    voxel: np.ndarray,
    origin: np.ndarray,
) -> float:
    idx = np.rint((points - origin) / voxel).astype(int)
    valid = np.all((idx >= 0) & (idx < np.array(data.shape[::-1])), axis=1)
    if valid.mean() < 0.99:
        return -1e9
    return float(data[idx[:, 2], idx[:, 1], idx[:, 0]].mean())


def write_transformed_pdb(
    source: Path,
    destination: Path,
    center: np.ndarray,
    target: np.ndarray,
    angle: float,
    reverse_z: bool,
) -> None:
    output = []
    for line in source.read_text().splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")):
            xyz = np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])]])
            x, y, z = transform(xyz, center, target, angle, reverse_z)[0]
            line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
        output.append(line)
    destination.write_text("".join(output))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("map", type=Path)
    parser.add_argument("pdb", type=Path)
    parser.add_argument("--out-pdb", type=Path, default=Path("analysis/mature_map/8FM9_coarse_fit.pdb"))
    parser.add_argument("--out-json", type=Path, default=Path("analysis/mature_map/8FM9_coarse_fit.json"))
    parser.add_argument(
        "--z-orientation",
        choices=("auto", "normal", "reversed"),
        default="auto",
        help="Constrain model Z orientation; reversed puts the 8FM9 floor toward increasing map Z.",
    )
    args = parser.parse_args()

    shape, voxel, origin, nsymbt = mrc_header(args.map)
    nx, ny, nz = shape
    data = np.memmap(args.map, dtype="<f4", mode="r", offset=1024 + nsymbt, shape=(nz, ny, nx))
    coords, groups = pdb_ca(args.pdb)
    center = coords.mean(axis=0)

    # The map is already centered in X/Y. Search a broad Z range and one
    # crystallographic asymmetric interval (30 degrees) around the C12 axis.
    map_xy = origin[:2] + voxel[:2] * np.array([(nx - 1) / 2, (ny - 1) / 2])
    coarse = []
    for reverse_z in (False, True):
        for angle in np.arange(0.0, 30.0, 1.0):
            for target_z in np.arange(250.0, 501.0, 2.5):
                target = np.array([map_xy[0], map_xy[1], target_z])
                value = score_points(transform(coords, center, target, angle, reverse_z), data, voxel, origin)
                coarse.append((value, reverse_z, float(angle), target.copy()))
    coarse.sort(key=lambda item: item[0], reverse=True)

    refined = []
    # Refine independent orientation hypotheses to avoid losing the second
    # possible membrane-facing direction too early.
    for reverse_z in (False, True):
        seed = next(item for item in coarse if item[1] == reverse_z)
        _, _, angle0, target0 = seed
        for angle in np.arange(angle0 - 2.0, angle0 + 2.01, 0.5):
            for dx in np.arange(-7.5, 7.51, 2.5):
                for dy in np.arange(-7.5, 7.51, 2.5):
                    for dz in np.arange(-7.5, 7.51, 2.5):
                        target = target0 + np.array([dx, dy, dz])
                        value = score_points(
                            transform(coords, center, target, angle, reverse_z), data, voxel, origin
                        )
                        refined.append((value, reverse_z, float(angle % 360), target.copy()))
    refined.sort(key=lambda item: item[0], reverse=True)
    orientation_best = {
        "normal": next(item for item in refined if not item[1]),
        "reversed": next(item for item in refined if item[1]),
    }
    if args.z_orientation == "normal":
        chosen = orientation_best["normal"]
    elif args.z_orientation == "reversed":
        chosen = orientation_best["reversed"]
    else:
        chosen = refined[0]
    best_score, reverse_z, angle, target = chosen
    fitted = transform(coords, center, target, angle, reverse_z)
    sampled = np.rint((fitted - origin) / voxel).astype(int)
    values = data[sampled[:, 2], sampled[:, 1], sampled[:, 0]]

    report = {
        "warning": "Provisional coarse rigid-body fit to a sharpened full map; not cross-validated.",
        "map": str(args.map),
        "source_pdb": str(args.pdb),
        "ca_atoms": int(len(coords)),
        "source_ca_center_xyz_angstrom": center.tolist(),
        "target_ca_center_xyz_angstrom": target.tolist(),
        "rotation_about_z_degrees": angle,
        "model_z_reversed": reverse_z,
        "mean_density_at_ca": best_score,
        "mean_density_at_floor_ca": float(values[groups == 0].mean()),
        "mean_density_at_pol_ca": float(values[groups == 1].mean()),
        "fraction_ca_above_map_mean_plus_2sd": float(
            np.mean(values > (float(data.mean()) + 2 * float(data.std())))
        ),
        "best_by_z_orientation": {
            name: {
                "score": float(item[0]),
                "reverse_z": bool(item[1]),
                "angle_deg": float(item[2]),
                "target_center_xyz_angstrom": item[3].tolist(),
            }
            for name, item in orientation_best.items()
        },
        "top_orientation_hypotheses": [
            {
                "score": float(item[0]),
                "reverse_z": bool(item[1]),
                "angle_deg": float(item[2]),
                "target_center_xyz_angstrom": item[3].tolist(),
            }
            for item in refined[:10]
        ],
    }
    args.out_pdb.parent.mkdir(parents=True, exist_ok=True)
    write_transformed_pdb(args.pdb, args.out_pdb, center, target, angle, reverse_z)
    args.out_json.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
