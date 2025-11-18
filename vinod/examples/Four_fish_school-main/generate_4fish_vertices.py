#!/usr/bin/env python3
"""
4-Fish Vertex Generator for IBAMR
Configuration: dx = 2.0L, dy = 0.4L

Formation:
   FISH-1: (x=0.0,  y=-0.2)  bottom-left   (Box 0)
   FISH-2: (x=2.0,  y=-0.2)  bottom-right  (Box 1)
   FISH-3: (x=0.0,  y=+0.2)  top-left      (Box 2)
   FISH-4: (x=2.0,  y=+0.2)  top-right     (Box 3)
"""

import os, sys

def transform_vertex_file(input_file, output_file, x_shift, y_shift):
    """Shift and recenter vertex coordinates."""
    with open(input_file, 'r') as f:
        n_vertices = int(f.readline().strip())
        coords = [list(map(float, line.split())) for line in f]
    y_mean = sum(y for _, y in coords) / len(coords)

    with open(output_file, 'w') as f:
        f.write(f"{n_vertices}\n")
        for x_old, y_old in coords:
            f.write(f"{x_old + x_shift:.6f}\t{y_old - y_mean + y_shift:.6f}\n")

def main():
    print("=" * 55)
    print("4-FISH VERTEX GENERATOR  (dx=2.0L, dy=0.4L)")
    print("=" * 55)

    base = "eel2d.vertex"
    backup = "eel2d_original.vertex"
    if not os.path.exists(base):
        print("[ERROR] Base file 'eel2d.vertex' not found.")
        sys.exit(1)
    if not os.path.exists(backup):
        os.system(f"cp {base} {backup}")
        print("[INFO] Created backup:", backup)

    fish_config = [
        ("eel2d_1.vertex", 0.0, -0.2),  # bottom-left
        ("eel2d_2.vertex", 2.0, -0.2),  # bottom-right
        ("eel2d_3.vertex", 0.0,  0.2),  # top-left
        ("eel2d_4.vertex", 2.0,  0.2),  # top-right
    ]

    for fname, dx, dy in fish_config:
        transform_vertex_file(backup, fname, dx, dy)
        print(f"[OK] {fname}  (x={dx}, y={dy})")

    print("\nAll vertex files generated successfully!")
    print("- eel2d_1.vertex → bottom-left  (Box 0)")
    print("- eel2d_2.vertex → bottom-right (Box 1)")
    print("- eel2d_3.vertex → top-left     (Box 2)")
    print("- eel2d_4.vertex → top-right    (Box 3)")
    print("=" * 55)

if __name__ == "__main__":
    main()
