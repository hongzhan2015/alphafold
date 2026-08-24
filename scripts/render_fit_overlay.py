#!/usr/bin/env python3
"""Render top and side projections of a fitted PDB over an MRC map."""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("map", type=Path)
    parser.add_argument("pdb", type=Path)
    parser.add_argument("--out", type=Path, default=Path("analysis/mature_map/8FM9_coarse_fit_overlay.png"))
    args = parser.parse_args()

    raw = args.map.read_bytes()[:1024]
    ints, floats = struct.unpack("<256i", raw), struct.unpack("<256f", raw)
    nx, ny, nz = ints[:3]
    voxel = np.array(floats[10:13]) / np.array(ints[7:10])
    origin = np.array(floats[49:52])
    data = np.memmap(args.map, dtype="<f4", mode="r", offset=1024 + ints[23], shape=(nz, ny, nx))

    coords, group = [], []
    for line in args.pdb.read_text().splitlines():
        if line.startswith("ATOM  ") and line[12:16].strip() == "CA":
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            group.append(0 if line[21] in "ABCDEFGHIJKL" else 1)
    indices = (np.asarray(coords) - origin) / voxel
    group = np.asarray(group)

    x0, x1, y0, y1, z0, z1 = 76, 212, 76, 212, 92, 196
    views = [
        (np.asarray(data[z0:z1, y0:y1, x0:x1].max(axis=0)), indices[:, [0, 1]], (x0, y0), "Top (XY)"),
        (np.asarray(data[z0:z1, y0:y1, x0:x1].max(axis=1)), indices[:, [0, 2]], (x0, z0), "Side (XZ)"),
        (np.asarray(data[z0:z1, y0:y1, x0:x1].max(axis=2)), indices[:, [1, 2]], (y0, z0), "Side (YZ)"),
    ]
    panels = []
    for array, points, offset, label in views:
        lo, hi = np.percentile(array, [5, 99.5])
        pixels = np.asarray(np.clip((array - lo) / (hi - lo), 0, 1) * 255, dtype=np.uint8)
        image = Image.fromarray(pixels, mode="L").convert("RGB").resize(
            (array.shape[1] * 3, array.shape[0] * 3), Image.Resampling.BILINEAR
        )
        draw = ImageDraw.Draw(image)
        for g, color in ((0, (0, 225, 255)), (1, (255, 190, 0))):
            q = points[group == g] - np.asarray(offset)
            for px, py in q:
                if 0 <= px < array.shape[1] and 0 <= py < array.shape[0]:
                    px, py = px * 3, py * 3
                    draw.ellipse((px - 1.2, py - 1.2, px + 1.2, py + 1.2), fill=color)
        canvas = Image.new("RGB", (image.width, image.height + 30), "white")
        canvas.paste(image, (0, 30))
        ImageDraw.Draw(canvas).text((8, 8), f"{label}: floor cyan; Pol gold", fill="black")
        panels.append(canvas)
    sheet = Image.new("RGB", (sum(p.width for p in panels), max(p.height for p in panels)), "white")
    x = 0
    for panel in panels:
        sheet.paste(panel, (x, 0))
        x += panel.width
    args.out.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(args.out)
    print(args.out)


if __name__ == "__main__":
    main()

