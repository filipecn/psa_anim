# this script generates the directory tree for a DSL simulation
# the files follow the same structure of faSavageHutterFOAM 1.0

import os
import argparse
import shutil
from pathlib import Path

def write_transport_properties(path, args, template_file):
    with open(template_file, "r") as f:
        file_content = f.read()
    file_content = file_content.replace("<<<EB>>>", str(args.eb));

    print("[init_dsl] writing ...", path / "transportProperties")
    with open(path / "transportProperties", "w") as f:
        f.write(file_content)

def write_control_dict_file(path, delta_t, write_time, end_time):
    control_file_path = path / "controlDict"
    with open(control_file_path, "r") as f:
        file_content = f.read()
    file_content = file_content.replace("<<<END_TIME>>>", "%.8f" % end_time) \
        .replace("<<<DELTA_T>>>", "%.8f" % delta_t) \
        .replace("<<<WRITE_INTERVAL>>>", "%.8f" % write_time)

    print("[init_dsl] writing ...", control_file_path)
    with open(control_file_path, "w") as f:
        f.write(file_content)


def write_all_run(path, mesh_generator, template_file, dem):
    with open(template_file, "r") as f:
        file_content = f.read()

    applications = ""
    if dem:
        applications += "runApplication gridToSTL\n"

    applications += "runApplication " + mesh_generator + "\n"

    file_content = file_content.replace("<<<APPLICATIONS>>>", str(applications))

    print("[init_dsl] writing ...", path / os.path.basename(template_file))
    with open(path / os.path.basename(template_file), "w") as f:
        f.write(file_content)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", type=Path, help="sim path")
    parser.add_argument("--template-dir", type=Path, help="template files dir")
    parser.add_argument("--dt", type=float, help="time step", default=0.0001)
    parser.add_argument("--wt", type=float, help="write time", default=0.2)
    parser.add_argument("--et", type=float, help="end time", default=90)
    parser.add_argument("--px", type=float, help="x domains", default=6)
    parser.add_argument("--py", type=float, help="y domains", default=1)
    parser.add_argument("--eb", type=float, help="erosion coeff", default=50)
    parser.add_argument("--snow-density", type=float, help="snow cloud density kg/m3", default=20)
    parser.add_argument("--mesh-generator", type=str, help="blockMesh|pMesh|cartesianMesh|makeFaMesh|slopeMesh",
                        default="cartesianMesh")
    parser.add_argument("--dem", action='store_true')
    args = parser.parse_args()

    # check templates
    if args.template_dir is None or not os.path.isdir(args.template_dir):
        exit(-1)

    # clean output directory
    if args.o.exists() and os.listdir(args.o):
        print("[init_dsl] Output folder not empty")
        exit(1)
    elif not args.o.exists():
        os.makedirs(args.o)

    constant_dir = args.o / "constant"
    system_dir = args.o / "system"

    # copy template files into new directory
    shutil.copytree(args.template_dir, args.o, symlinks=False, ignore=None, ignore_dangling_symlinks=False,
                    dirs_exist_ok=True)

    # file template parameters
    write_transport_properties(constant_dir, args, args.template_dir / "constant/transportProperties")
    write_control_dict_file(system_dir, args.dt, args.wt, args.et)
    write_all_run(args.o, args.mesh_generator, args.template_dir / "Allrun", args.dem)
    write_all_run(args.o, args.mesh_generator, args.template_dir / "Allrun-parallel", args.dem)
