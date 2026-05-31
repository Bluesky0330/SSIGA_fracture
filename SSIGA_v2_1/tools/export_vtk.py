#!/usr/bin/env python3
import os

ROOT = os.path.dirname(os.path.dirname(__file__))
INPUT_VERT = os.path.join(ROOT, 'analysis', 'input', 'SCCHP_v10_coarse', 'input_vertex.txt')
INPUT_CTRL = os.path.join(ROOT, 'analysis', 'output_SCCHP_v10_coarse', 'second_inverse_new_control_points.txt')
OUT_DIR = os.path.join(ROOT, 'analysis')


def read_vertices(path):
    pts = []
    elem_ids = []
    point_ids = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.lower().startswith('element_id'):
                continue
            parts = line.split()
            if len(parts) < 5:
                # try tab split
                parts = line.split('\t')
            if len(parts) < 5:
                continue
            e = int(float(parts[0]))
            pid = int(float(parts[1]))
            x = float(parts[2])
            y = float(parts[3])
            z = float(parts[4])
            pts.append((x,y,z))
            elem_ids.append(e)
            point_ids.append(pid)
    return pts, elem_ids, point_ids


def read_control_points(path):
    pts = []
    weights = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            # format: id x y [z] w
            if len(parts) == 5:
                _id, x, y, z, w = parts
                x, y, z = float(x), float(y), float(z)
                w = float(w)
            elif len(parts) == 4:
                _id, x, y, w = parts
                x, y = float(x), float(y)
                z = 0.0
                w = float(w)
            else:
                # more columns -> assume id x y z w ...
                _id = parts[0]
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                w = float(parts[4]) if len(parts) > 4 else 1.0
            pts.append((x,y,z))
            weights.append(w)
    return pts, weights


def write_vtk_points(path, pts, arrays=None, array_names=None):
    n = len(pts)
    with open(path, 'w', encoding='utf-8') as f:
        f.write('# vtk DataFile Version 3.0\n')
        f.write(os.path.basename(path) + '\n')
        f.write('ASCII\n')
        f.write('DATASET POLYDATA\n')
        f.write(f'POINTS {n} float\n')
        for x,y,z in pts:
            f.write(f'{x:.15e} {y:.15e} {z:.15e}\n')
        if arrays and array_names:
            f.write(f'\nPOINT_DATA {n}\n')
            for name, arr in zip(array_names, arrays):
                # determine type
                if all(isinstance(v, int) for v in arr):
                    f.write(f'SCALARS {name} int 1\nLOOKUP_TABLE default\n')
                    for v in arr:
                        f.write(f'{v}\n')
                else:
                    f.write(f'SCALARS {name} float 1\nLOOKUP_TABLE default\n')
                    for v in arr:
                        f.write(f'{float(v):.15e}\n')


if __name__ == '__main__':
    print('Reading vertices from', INPUT_VERT)
    pts, elem_ids, point_ids = read_vertices(INPUT_VERT)
    print('Read', len(pts), 'points')
    out_vertices = os.path.join(OUT_DIR, 'input_vertex.vtk')
    write_vtk_points(out_vertices, pts, arrays=[elem_ids, point_ids], array_names=['element_id','point_id'])
    print('Wrote', out_vertices)

    if os.path.exists(INPUT_CTRL):
        print('Reading control points from', INPUT_CTRL)
        cpts, weights = read_control_points(INPUT_CTRL)
        out_ctrl = os.path.join(OUT_DIR, 'control_points.vtk')
        write_vtk_points(out_ctrl, cpts, arrays=[weights], array_names=['w'])
        print('Wrote', out_ctrl)
    else:
        print('Control points file not found:', INPUT_CTRL)
