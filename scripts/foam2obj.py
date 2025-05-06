# This script exports the openfoam geometry into an obj file

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

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        type=Path,
        help="simulation directory containing constant/polyMesh",
        required=True,
    )
    parser.add_argument("-o", type=Path, help="output directory")
    parser.add_argument("-p", type=str, help="patch name")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    verbose = args.verbose

    LOG("[foam2obj] Starting foam2obj --------------------------")

    # check paths
    if args.i is None or not os.path.isdir(args.i):
        ERR("[foam2obj] invalid input directory")
        ERR("[foam2obj]", args.i)
        exit(-1)

    output_path = args.o

    if args.p is not None:
        output_obj = output_path / (args.p + ".obj")

        # check if output file exists
        if os.path.isfile(output_obj):
            ERR("[foam2obj] output file already exists")
            exit(1)

        LOG("[foam2obj] input:  ", str(args.i))
        LOG("[foam2obj] output: ", str(output_obj))
        LOG("[foam2obj] patch:  ", args.p)

        # mesh data
        of_sim = psa_anim_py.OFSim(args.verbose)
        of_sim.setSimPath(str(args.i))

        # patch = of_sim.DSL(args.p)

        of_sim.exportOBJ(args.p, str(output_obj))

    LOG("[foam2obj] ------------------------ complete")
