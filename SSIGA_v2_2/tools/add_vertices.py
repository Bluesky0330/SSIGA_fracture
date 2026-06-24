import sys
from pathlib import Path

in_path = Path('analysis/input_vertex.vtk')
out_path = Path('analysis/input_vertex_with_vertices.vtk')
if not in_path.exists():
    print('input file not found:', in_path)
    sys.exit(1)

with in_path.open('r', encoding='utf-8') as fin, out_path.open('w', encoding='utf-8') as fout:
    # copy header until POINTS line
    for line in fin:
        fout.write(line)
        if line.strip().upper().startswith('POINTS '):
            parts = line.split()
            if len(parts) < 2:
                print('unable to parse POINTS header')
                sys.exit(1)
            n_points = int(parts[1])
            # copy next n_points coordinate lines
            for i in range(n_points):
                l = fin.readline()
                if not l:
                    print('unexpected EOF while reading points')
                    sys.exit(1)
                fout.write(l)
            # write VERTICES section
            total_indices = n_points + n_points  # n + sum(counts)=n + n*1
            fout.write(f'VERTICES {n_points} {total_indices}\n')
            # write indices in lines of reasonable length
            for i in range(n_points):
                fout.write(f'1 {i}\n')
            # now copy the rest of the file
            break
    # copy remaining
    for line in fin:
        fout.write(line)

print('wrote', out_path)
