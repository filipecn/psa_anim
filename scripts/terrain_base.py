#!/bin/bash
# this script takes a surface mesh and extruds downards to create a solid base terrain

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
    print("[terrain_base] using psa_anim_py from", psa_anim_py_path)

sys.path.insert(0, psa_anim_py_path)
import psa_anim_py

verbose = False


def LOG(*args):
    if verbose:
        print(args)


def ERR(*args):
    print(args)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=Path, help="mesh file", required=True)
    parser.add_argument("-o", type=Path, help="output obj", required=True)
    parser.add_argument("--base-height", type=float, help="base height", default=100)
    parser.add_argument(
        "--bottom-type", type=str, help="[xy | inclined | extrude]", default="xy"
    )
    parser.add_argument("--close-top", action="store_true", help="include input mesh")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    verbose = args.verbose

    LOG("[terrain_base] Starting terrain_base --------------------------")

    if not os.path.exists(args.i):
        ERR("[terrain_base] Invalid file path")
        exit(1)

    if os.path.exists(args.o):
        ERR("[terrain_base] Output file exists")
        exit(1)

    if args.bottom_type == "xy":
        bottom_type = 0
    elif args.bottom_type == "inclined":
        bottom_type = 1
    elif args.bottom_type == "extrude":
        bottom_type = 2
    else:
        ERR("[terrain_base] Invalid bottom type option")
        exit(1)

    LOG("[terrain_base] input:       ", str(args.i))
    LOG("[terrain_base] output:      ", str(args.o))
    LOG("[terrain_base] base height: ", args.base_height)
    LOG("[terrain_base] bottom type: ", args.bottom_type)
    LOG("[terrain_base] close top:   ", args.close_top)

    psa_anim_py.terrain_base(
        str(args.i),
        str(args.o),
        args.base_height,
        bottom_type,
        args.close_top,
        args.verbose,
    )

    LOG("[terrain_base] ------------------------ complete")
