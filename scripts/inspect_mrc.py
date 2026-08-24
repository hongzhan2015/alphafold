#!/usr/bin/env python3
"""Inspect an MRC/CCP4 map using only NumPy and Pillow."""

from __future__ import annotations

import argparse
import json
import struct
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


MODE_DTYPES = {0: "i1", 1: "<i2", 2: "<f4", 6: "<u2"}


def read_header(path: Path) -> dict:
    raw = path.read_bytes()[:1024]
    ints = struct.unpack("<256i", raw)
    floats = struct.unpack("<256f", raw)
    nx, ny, nz, mode = ints[:4]
    if mode not in MODE_DTYPES:
        raise ValueError(f"Unsupported MRC mode {mode}")
    mx, my, mz = ints[7:10]
    cell = floats[10:13]
    voxel = [cell[0] / mx, cell[1] / my, cell[2] / mz]
    return {
        "shape_xyz": [nx, ny, nz],
        "shape_zyx": [nz, ny, nx],
        "mode": mode,
        "dtype": MODE_DTYPES[mode],
        "start_xyz": list(ints[4:7]),
        "grid_xyz": [mx, my, mz],
        "cell_angstrom": list(cell),
        "voxel_size_angstrom_xyz": voxel,
        "axis_order_mapc_mapr_maps": list(ints[16:19]),
        "header_density_min_max_mean": list(floats[19:22]),
        "ispg": ints[22],
        "nsymbt": ints[23],
        "origin_xyz_angstrom": list(floats[49:52]),
        "header_rms": floats[54],
        "map_signature": raw[208:212].decode("ascii", "replace"),
    }


def as_uint8(array: np.ndarray, low: float, high: float) -> np.ndarray:
    scaled = np.clip((array.astype(np.float32) - low) / (high - low), 0, 1)
    return np.asarray(np.rint(scaled * 255), dtype=np.uint8)


def labeled_panel(array: np.ndarray, label: str, low: float, high: float) -> Image.Image:
    image = Image.fromarray(as_uint8(array, low, high), mode="L").convert("RGB")
    image = image.resize((432, 432), Image.Resampling.BILINEAR)
    canvas = Image.new("RGB", (432, 462), "white")
    canvas.paste(image, (0, 30))
    ImageDraw.Draw(canvas).text((8, 8), label, fill="black")
    return canvas


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("map", type=Path)
    parser.add_argument("--out-dir", type=Path, default=Path("analysis/mature_map"))
    args = parser.parse_args()

    header = read_header(args.map)
    nx, ny, nz = header["shape_xyz"]
    offset = 1024 + header["nsymbt"]
    data = np.memmap(
        args.map,
        dtype=np.dtype(header["dtype"]),
        mode="r",
        offset=offset,
        shape=(nz, ny, nx),
        order="C",
    )
    finite = np.asarray(data[np.isfinite(data)], dtype=np.float32)
    nonzero = finite[finite != 0]
    source = nonzero if nonzero.size else finite
    percentiles = [0.1, 0.5, 1, 5, 25, 50, 75, 95, 99, 99.5, 99.9]
    stats = {
        "file": str(args.map),
        "file_size_bytes": args.map.stat().st_size,
        "finite_voxels": int(finite.size),
        "nonzero_voxels": int(nonzero.size),
        "zero_fraction": float(1 - nonzero.size / finite.size),
        "actual_min": float(finite.min()),
        "actual_max": float(finite.max()),
        "actual_mean": float(finite.mean()),
        "actual_std": float(finite.std()),
        "nonzero_percentiles": {
            str(p): float(v) for p, v in zip(percentiles, np.percentile(source, percentiles))
        },
    }

    threshold = stats["actual_mean"] + 2.0 * stats["actual_std"]
    mask = np.asarray(data > threshold)
    coords = np.argwhere(mask)
    if coords.size:
        weights = np.asarray(data[mask], dtype=np.float64) - threshold
        center_zyx = np.average(coords, axis=0, weights=weights)
        bbox_min = coords.min(axis=0)
        bbox_max = coords.max(axis=0)
    else:
        center_zyx = np.array([nz / 2, ny / 2, nx / 2])
        bbox_min = np.zeros(3, dtype=int)
        bbox_max = np.array([nz - 1, ny - 1, nx - 1])
    center_indices = np.rint(center_zyx).astype(int)
    stats["two_sigma_threshold"] = float(threshold)
    stats["two_sigma_voxels"] = int(mask.sum())
    stats["positive_density_center_zyx_voxels"] = center_zyx.tolist()
    stats["positive_density_center_xyz_angstrom"] = [
        float(center_zyx[2] * header["voxel_size_angstrom_xyz"][0]),
        float(center_zyx[1] * header["voxel_size_angstrom_xyz"][1]),
        float(center_zyx[0] * header["voxel_size_angstrom_xyz"][2]),
    ]
    stats["two_sigma_bbox_zyx_inclusive"] = [bbox_min.tolist(), bbox_max.tolist()]

    z, y, x = center_indices
    slab = 12
    arrays = [
        (np.asarray(data[z]), f"XY slice at z={z}"),
        (np.asarray(data[:, y, :]), f"XZ slice at y={y}"),
        (np.asarray(data[:, :, x]), f"YZ slice at x={x}"),
        (
            np.asarray(data[max(0, z-slab):min(nz, z+slab+1)].max(axis=0)),
            f"XY max projection, z={z-slab}..{z+slab}",
        ),
        (
            np.asarray(data[:, max(0, y-slab):min(ny, y+slab+1), :].max(axis=1)),
            f"XZ max projection, y={y-slab}..{y+slab}",
        ),
        (
            np.asarray(data[:, :, max(0, x-slab):min(nx, x+slab+1)].max(axis=2)),
            f"YZ max projection, x={x-slab}..{x+slab}",
        ),
    ]
    low = float(np.percentile(source, 5))
    high = float(np.percentile(source, 99.7))
    panels = [labeled_panel(a, label, low, high) for a, label in arrays]
    sheet = Image.new("RGB", (432 * 3, 462 * 2), "white")
    for index, panel in enumerate(panels):
        sheet.paste(panel, ((index % 3) * 432, (index // 3) * 462))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "map_report.json").write_text(
        json.dumps({"header": header, "statistics": stats}, indent=2) + "\n"
    )
    sheet.save(args.out_dir / "map_views.png")
    print(json.dumps({"header": header, "statistics": stats}, indent=2))
    print(args.out_dir / "map_views.png")


if __name__ == "__main__":
    main()

