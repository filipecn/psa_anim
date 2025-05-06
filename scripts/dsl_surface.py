# This script generates a polygonal mesh of a the DSL simulation output
# A single mesh is generated for each simulation time step and uses the height field
# The terrain surface is extruded in the normal direction based on the height information

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
    print("[dsl_surface] using psa_anim_py from", psa_anim_py_path)
sys.path.insert(0, psa_anim_py_path)
import psa_anim_py

verbose = False


def LOG(*args):
    if verbose:
        print(args)


def ERR(*args):
    print(args)


if __name__ == "__main__":
    LOG("[dsl_surface] ---------------------------------- start")

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        type=Path,
        help="dsl openfoam directory (containing constant/polymesh)",
        required=True,
    )
    parser.add_argument("-o", type=Path, help="output directory", required=True)
    parser.add_argument("-p", type=str, help="terrain patch name", required=True)
    parser.add_argument("-s", type=Path, help="dsl simulation data directory")
    parser.add_argument("-F", action="store_true", help="force write")
    parser.add_argument(
        "--only-missing-frames", action="store_true", help="write only missing frames"
    )
    parser.add_argument(
        "--with-height-map", action="store_true", help="write height colors"
    )
    parser.add_argument(
        "--with-distance-map", action="store_true", help="write front distance colors"
    )
    parser.add_argument("-f", type=int, help="single frame number")
    parser.add_argument(
        "--renumber-frames",
        action="store_true",
        help="re-number frames to a contiguous sequence",
        default=False,
    )
    parser.add_argument(
        "--max-height", type=float, help="max allowed height (for coloring)", default=4
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    verbose = args.verbose

    # check paths
    if args.i is None or not os.path.isdir(args.i):
        ERR("[dsl_surface] Input directory could not be found!")
        exit()

    if args.f is not None and not os.path.isdir(args.f):
        ERR("[dsl_surface] Input frames directory could not be found!")
        exit()

    # parameters
    input_dir = args.i
    single_frame = args.f or -1
    frames_dir = args.s or input_dir

    output_path = args.o
    if not output_path.exists():
        os.makedirs(output_path)

    LOG("[dsl_surface] input:        ", str(input_dir))
    LOG("[dsl_surface] output:       ", str(output_path))
    LOG("[dsl_surface] frames dir:   ", str(frames_dir))
    LOG("[dsl_surface] single frame: ", str(single_frame))

    # simulation data
    of_sim = psa_anim_py.OFSim(args.verbose)
    of_sim.setSimPath(str(input_dir))
    of_sim.setFramesPath(str(frames_dir))

    # get DSL
    patch = of_sim.DSL(args.p, args.verbose)

    # gen surface
    patch.extractSurface(
        str(output_path),
        args.max_height,
        args.with_height_map,
        args.with_distance_map,
        args.renumber_frames,
        args.only_missing_frames,
        args.F,
        single_frame,
    )

    LOG("[dsl_surface] ---------------------------------- end")
