# This script projects a quad-mesh grid onto a input mesh

import sys
import os
import argparse
from pathlib import Path

psa_anim_py_path = os.getenv("PSA_ANIM_PY_PATH")
if psa_anim_py_path is None:
    print(
        "Please set PSA_ANIM_PY_PATH environment variable! (run 'source <build>/psa_anim_variables.sh')"
    )
    exit(-1)
else:
    print("[mesh2quad] using psa_anim_py from", psa_anim_py_path)

sys.path.insert(0, psa_anim_py_path)
import psa_anim_py


def compute_axis(stl_file):
    x_bounds = [0, 0]
    y_bounds = [0, 0]
    first = True
    print(stl_file)
    with open(stl_file, "r") as stl:
        lines = stl.readlines()
        for line in lines:
            words = line.split()
            if words[0] == "vertex":
                x = float(words[1])
                y = float(words[2])
                if first:
                    first = False
                    print(x, y, line)
                x_bounds[0] = min([x, x_bounds[0]])
                x_bounds[1] = max([x, x_bounds[1]])
                y_bounds[0] = min([y, y_bounds[0]])
                y_bounds[1] = max([y, y_bounds[1]])

    return [
        x_bounds[0],
        x_bounds[1],
        (y_bounds[0] + y_bounds[1]) * 0.5,
        y_bounds[1] - y_bounds[0],
    ]


if __name__ == "__main__":
    print("[mesh2quad] Starting mesh2quad --------------------------")

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=Path, help=".stl file", required=True)
    parser.add_argument("-o", type=Path, help="output obj", required=True)
    parser.add_argument(
        "--solid-name", type=str, help="surface patch name", default="terrain"
    )
    parser.add_argument("--axis-a", type=float, nargs="+", help="[x, y, z]")
    parser.add_argument("--axis-b", type=float, nargs="+", help="[x, y, z]")
    parser.add_argument("--width", type=float, help="grid lateral width")
    parser.add_argument("--cell-size", type=float, help="cell size", default=5)
    parser.add_argument(
        "--terrain-aligned", action="store_true", help="descend is x ligned"
    )
    parser.add_argument("--d2", action="store_true")
    args = parser.parse_args()

    if not os.path.exists(args.i):
        print("[mesh2quad] Invalid stl file path")
        print("[mesh2quad][stl file] " + str(args.i))
        exit(1)

    if os.path.exists(args.o):
        print("[mesh2quad] Output file exists")

    print("[mesh2quad] input:       ", args.i)
    print("[mesh2quad] output:      ", args.o)

    if args.terrain_aligned:
        bounds = compute_axis(args.i)
        axis_a = [bounds[1], bounds[2]]
        axis_b = [bounds[0], bounds[2]]
        width = bounds[3]
    else:
        axis_a = args.axis_a
        axis_b = args.axis_b
        width = args.width

    print("[mesh2quad] axis points: ", axis_a, axis_b)
    print("[mesh2quad] cell size:   ", args.cell_size)
    print("[mesh2quad] width:       ", width, flush=True)

    psa_anim_py.convert_to_quad_mesh(
        str(args.i), str(args.o), axis_a, axis_b, args.cell_size, width, args.d2
    )

    print("[mesh2quad] ------------------------ complete")
