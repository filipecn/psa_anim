# This script prepares and computes the input data for the PSL simulation
# The input data is computed from the DSL simulation output
# Once the DSL is fully simulated and its data is stored, the following
# fields are computed/extracted:
# - dsl_h            extracted
# - dsl_U            extracted
# - slope            computed
# - dsl_front_sdf        computed
#
# The fields listed above belong to the surface patch (terrain, bottom, etc...)
# and the output fields are also surface fields for the same patch.
# The output files must be stored under the path constant/boundaryData of the
# PSL simulation directory, so openfoam can access them properly.
#
# This script receives the following inputs:
# - DSL directory
# - PSL directory
# - DSL patch name (usually 'terrain')
# - PSL patch name (usually 'bottom')

import sys
import os
import argparse
from pathlib import Path

psa_anim_py_path = os.getenv("PSA_ANIM_PY_PATH")
if psa_anim_py_path is None:
    print(
        "Please set PSA_ANIM_PY_PATH environment variable!"
        "(run 'source <build>/psa_anim_variables.sh')"
    )
    exit(-1)
else:
    print("[dsl2psl] using psa_anim_py from", psa_anim_py_path)

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
    parser.add_argument("-i", type=Path, help="simulation directory")
    parser.add_argument("-o", type=Path, help="output directory")
    parser.add_argument("-p", type=str, help="input patch name")
    parser.add_argument("-op", type=str, help="output patch name")
    parser.add_argument("-f", type=Path, help="frames path")
    parser.add_argument("-j", action="store_true", help="jitter points")
    parser.add_argument("--empty-dsl", action="store_true", help="generate empty dsl")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    verbose = args.verbose

    # check paths
    if args.i is None or not os.path.isdir(args.i):
        ERR("[dsl2psl] input path not found")
        exit(1)

    output_path = args.o / args.op

    if output_path.exists() and os.listdir(output_path):
        ERR("[dsl2psl] output folder not empty")
        exit(1)
    elif not output_path.exists():
        os.makedirs(output_path)

    LOG("[dsl2psl] convert dsl data into psl boundary data")
    LOG("[dsl2psl] ---------------------------------------")
    LOG("[dsl2psl][simulation directory] " + str(args.i))
    LOG("[dsl2psl][output directory]     " + str(output_path))
    LOG("[dsl2psl][input patch name]     " + str(args.p))
    LOG("[dsl2psl][output patch name]    " + str(args.op))
    LOG("[dsl2psl][frames path]          " + str(args.f))
    LOG("[dsl2psl][jitter points]        " + str(args.j))
    LOG("[dsl2psl][generate empty dsl]   " + str(args.empty_dsl))

    # simulation data
    of_sim = psa_anim_py.OFSim(args.verbose)
    of_sim.setSimPath(str(args.i))

    if args.f is not None:
        of_sim.setFramesPath(str(args.f))

    # get DSL
    patch = of_sim.DSL(args.p, args.verbose)

    # produce PSL Input
    patch.preparePSLInput(str(output_path), args.j, args.empty_dsl)

    LOG("[dsl2psl] ------------------------ complete")
