"""Extract the fourteen black theory profiles from the published vector PDF.

Usage: python extract_fig3.py Khadir_JOSAA_19.pdf
Requires Poppler's pdftocairo and NumPy. Tests use the saved CSV offline.
This measures the published curves, independently of the simulation code.
"""

from pathlib import Path
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

import numpy as np


def outline_points(d):
    tokens = re.findall(r'[MLCZ]|-?\d+(?:\.\d+)?', d)
    points, i = [], 0
    while i < len(tokens):
        command = tokens[i]
        i += 1
        if command in ('M', 'L'):
            point = np.array(tokens[i:i + 2], dtype=float)
            i += 2
            points.append(point)
            if command == 'M':
                start = point
        elif command == 'C':
            p0 = points[-1]
            p1, p2, p3 = np.array(tokens[i:i + 6], dtype=float).reshape(3, 2)
            i += 6
            for t in np.linspace(0, 1, 9)[1:]:
                points.append((1 - t)**3 * p0 + 3 * (1 - t)**2 * t * p1 + 3 * (1 - t) * t**2 * p2 + t**3 * p3)
        elif command == 'Z':
            points.append(start)
    return np.asarray(points)


def extract(pdf):
    with tempfile.TemporaryDirectory() as directory:
        svg = Path(directory) / 'page.svg'
        subprocess.run(['pdftocairo', '-f', '4', '-l', '4', '-svg', str(pdf), str(svg)], check=True)
        paths = ET.parse(svg).getroot().findall('.//{http://www.w3.org/2000/svg}path')
    curves, gridlines = [], []
    for path in paths:
        if path.get('fill') != 'rgb(0%, 0%, 0%)':
            continue
        numbers = np.array(re.findall(r'-?\d+(?:\.\d+)?', path.get('d', '')), dtype=float)
        if not len(numbers) or len(numbers) % 2:
            continue
        x, y = numbers[::2], numbers[1::2]
        if x.min() < 100 or x.max() > 495 or y.min() < 430 or y.max() > 685:
            continue
        if len(path.get('d', '')) > 5000 and y.max() - y.min() > 1:
            curves.append(outline_points(path.get('d')))
        if 48 < x.max() - x.min() < 52 and y.max() - y.min() < 0.3:
            gridlines.append(((x.min() + x.max()) / 2, (y.min() + y.max()) / 2))
    assert len(curves) == 14
    curves.sort(key=lambda points: (points[:, 1].mean() > 600, points[:, 0].mean()))
    spacing = 2 * np.pi / (801 * (1.5 * 2 * np.pi / 0.590 / 100))
    x_object = np.linspace(-4.5, 4.5, 181)
    rows = []
    z = [-1.4, -0.98, -0.55, -0.02, 0.55, 0.94, 1.4]
    for index, points in enumerate(curves):
        x, y = points.T
        center = (x.min() + x.max()) / 2
        is_opd = index >= 7
        grid_y = sorted(gy for gx, gy in gridlines if abs(gx - center) < 2 and (gy > 600) == is_opd)
        levels = [10, 5, 0, -5, -10] if is_opd else [1.05, 1.0, 0.95, 0.90, 0.85, 0.80]
        assert len(grid_y) == len(levels)
        origin = 0.0 if is_opd else 1.0
        zero_y = grid_y[levels.index(origin)]
        offsets = np.asarray(grid_y) - zero_y
        slope = np.dot(offsets, np.asarray(levels) - origin) / np.dot(offsets, offsets)
        calibration = [slope, origin - slope * zero_y]
        # The line is a filled 0.25-pt stroke around 201 data samples.
        x_pdf = center + x_object / (100 * spacing) * (x.max() - x.min() - 0.25) / 2
        p0, p1 = points[:-1], points[1:]
        values = []
        for xp in x_pdf:
            crosses = (p0[:, 0] <= xp) != (p1[:, 0] <= xp)
            a, b = p0[crosses], p1[crosses]
            intersections = a[:, 1] + (xp - a[:, 0]) * (b[:, 1] - a[:, 1]) / (b[:, 0] - a[:, 0])
            values.append(np.polyval(calibration, (intersections.min() + intersections.max()) / 2))
        rows.extend((int(is_opd), z[index % 7], xp, value) for xp, value in zip(x_object, values))
    np.savetxt(Path(__file__).with_name('khadir_fig3_published.csv'), rows, delimiter=',',
               header='observable_0_T_1_OPDnm,z_um,x_object_um,value', fmt='%.9g')


if __name__ == '__main__':
    extract(Path(sys.argv[1]))
