# This script takes quad mesh from a .obj file and converts it into a OpenFOAM's block mesh dictionary

import sys
import os
import argparse
from pathlib import Path
import math

psa_anim_py_path = os.getenv("PSA_ANIM_PY_PATH")
if psa_anim_py_path is None:
    print(
        "Please set PSA_ANIM_PY_PATH environment variable! (run 'source <build>/psa_anim_variables.sh')"
    )
    exit(-1)
else:
    print("[quad2block_mesh_desc] using psa_anim_py from", psa_anim_py_path)

sys.path.insert(0, psa_anim_py_path)
import psa_anim_py


class Vec:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def normalized(self):
        s = self.x**2 + self.y**2
        d = math.sqrt(s)
        return Vec(self.x / d, self.y / d)

    def __sub__(self, other):
        return Vec(self.x - other.x, self.y - other.y)

    def left(self):
        return Vec(-self.y, self.x)

    def right(self):
        return Vec(self.y, -self.x)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=Path, help=".obj file", required=True)
    parser.add_argument(
        "-o", type=Path, help="output directory (system folder)", required=True
    )
    parser.add_argument("--grading", type=float, nargs="+", help="[x, y, z]")
    parser.add_argument("--res", type=int, nargs="+", help="[x, y, z]")
    parser.add_argument("--inclined-top", action="store_true")
    parser.add_argument("--height", type=float, default=50)
    parser.add_argument(
        "--axis-a", type=float, nargs="+", help="[x, y, z]", default=[0, 0]
    )
    parser.add_argument(
        "--axis-b", type=float, nargs="+", help="[x, y, z]", default=[1, 0]
    )
    parser.add_argument("--d2", action="store_true")
    parser.add_argument("--domain-type", type=str, help="openBox", default="openBox")
    args = parser.parse_args()

    # check paths
    if args.i is None or not os.path.isfile(args.i):
        print("[quad2block_mesh_desc] Invalid obj file")
        exit(1)

    if args.o.exists():
        print("[quad2block_mesh_desc] Output file not empty")
        exit(1)

    print("[quad2block_mesh_desc] Starting quad2block_mesh_desc -------------")
    print("[quad2block_mesh_desc] input:         ", str(args.i))
    print("[quad2block_mesh_desc] outout:        ", str(args.o))
    print("[quad2block_mesh_desc] grading:       ", args.grading)
    print("[quad2block_mesh_desc] res:           ", args.res)
    print("[quad2block_mesh_desc] inclined-top:  ", args.inclined_top)
    print("[quad2block_mesh_desc] height:        ", args.height)
    print("[quad2block_mesh_desc] axis-a:        ", args.axis_a)
    print("[quad2block_mesh_desc] axis-b:        ", args.axis_b)
    print("[quad2block_mesh_desc] domain-type:   ", args.domain_type)
    print("[quad2block_mesh_desc] d2:            ", args.d2)

    output_path = args.o

    a = Vec(args.axis_a[0], args.axis_a[1])
    b = Vec(args.axis_b[0], args.axis_b[1])
    x_ = (b - a).normalized()
    y_ = x_.right()

    bm = psa_anim_py.BlockMeshDesc()

    bm.addPatchDirection("top", [0, 0, 1])
    bm.addPatchDirection("terrain", [0, 0, -1])

    bm.addPatchDirection("outlet", [x_.x, x_.y, 0])
    bm.addPatchDirection("inlet", [-x_.x, -x_.y, 0])

    bm.addPatchDirection("walls", [y_.x, y_.y, 0])
    bm.addPatchDirection("walls", [-y_.x, -y_.y, 0])

    if args.domain_type == "tunnel":
        bm.setBoundaryType("top", "wall")
        bm.setBoundaryType("terrain", "wall")
        bm.setBoundaryType("inlet", "patch")
        bm.setBoundaryType("outlet", "patch")
        bm.setBoundaryType("walls", "wall")

    if args.d2:
        bm.setBoundaryType("top", "patch")
        bm.setBoundaryType("terrain", "wall")
        bm.setBoundaryType("inlet", "patch")
        bm.setBoundaryType("outlet", "patch")
        bm.setBoundaryType("walls", "empty")

    if args.res is not None:
        bm.setBlockResolution(args.res)

    if args.grading is not None:
        bm.setBlockGrading(args.grading)

    bm.setSlopeDirection(
        [args.axis_b[0] - args.axis_a[0], args.axis_b[1] - args.axis_a[1]]
    )
    bm.setSlopePoint(args.axis_a)
    bm.setHeight(args.height)
    bm.setTopIsInclined(args.inclined_top)
    bm.loadOBJ(str(args.i))
    bm.save(str(args.o))

    print("[quad2block_mesh_desc] ------------------------ complete")
