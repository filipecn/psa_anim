# this script creates the directory tree for a PSL simulation
# the following directories and files are created
# sim path {
#   constant {
#       - transportProperties
#   }
#   system {
#       - controlDict
#       - decomposeParDict
#       - fvSchemes
#       - fvSolution
#       - pslControl
#   }
# }

import os
import argparse
import shutil
from pathlib import Path

foam_header = (
    "/*--------------------------------*- C++ -*----------------------------------*\n"
    "| =========                 |                                                 |\n"
    "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
    "|  \\    /   O peration     | Version:  v2012                                 |\n"
    "|   \\  /    A nd           | Website:  www.openfoam.com                      |\n"
    "|    \\/     M anipulation  |                                                 |\n"
    "\\*---------------------------------------------------------------------------*/\n"
)
foam_footer = (
    "// ************************************************************************* //\n"
)

foam_file_dict_template = (
    "FoamFile\n{\n"
    "\tversion     2.0;\n"
    "\tformat      ascii;\n"
    "\tclass       dictionary;\n"
    "%s\n}\n"
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
)


def write_control_dict_file(
    path, delta_t, write_time, end_time, write_format, max_co, template_file
):
    with open(template_file, "r") as f:
        file_content = f.read()
    file_content = (
        file_content.replace("<<<APPLICATION>>>", "pslFoam")
        .replace("<<<END_TIME>>>", "%.8f" % end_time)
        .replace("<<<DELTA_T>>>", "%.8f" % delta_t)
        .replace("<<<WRITE_FORMAT>>>", write_format)
        .replace("<<<WRITE_INTERVAL>>>", "%.8f" % write_time)
        .replace("<<<MAX_CO>>>", "%.8f" % max_co)
    )

    print("[init_psl] writing ...", path / "controlDict")
    with open(path / "controlDict", "w") as f:
        f.write(file_content)


def write_transport_properties(
    path, snow_density, air_density, template_file, dab, snow_viscosity
):
    with open(template_file, "r") as f:
        file_content = f.read()
    file_content = (
        file_content.replace("<<<SNOW_DENSITY>>>", str(snow_density))
        .replace("<<<SNOW_VISCOSITY>>>", str(snow_viscosity))
        .replace("<<<AIR_DENSITY>>>", str(air_density))
        .replace("<<<DAB>>>", str(dab))
    )

    print("[init_psl] writing ...", path / "transportProperties")
    with open(path / "transportProperties", "w") as f:
        f.write(file_content)


def write_turbulence_properties(path, use_turbulence, template_file):
    with open(template_file, "r") as f:
        file_content = f.read()

    if use_turbulence:
        file_content = file_content.replace("<<<SIMULATION_TYPE>>>", "RAS").replace(
            "<<<CONFIG>>>",
            "RAS\n"
            "{\n"
            "  // Mandatory entries\n"
            "  RASModel        kEpsilon;\n"
            "  // Optional entries\n"
            "  turbulence      on;\n"
            "  printCoeffs     on;\n"
            "  // Optional model coefficients\n"
            "  Cmu             0.09;\n"
            "  C1              1.44;\n"
            "  C2              1.92;\n"
            "  C3              0.0;\n"
            "  sigmak          1.0;\n"
            "  sigmaEps        1.3;\n"
            "}\n",
        )
    else:
        file_content = file_content.replace("<<<SIMULATION_TYPE>>>", "laminar").replace(
            "<<<CONFIG>>>", ""
        )

    print("[init_psl] writing ...", path / "turbulenceProperties")
    with open(path / "turbulenceProperties", "w") as f:
        f.write(file_content)


def write_fv_schemes(path, use_turbulence, template_file):
    div_schemes = {
        "div(rhoPhi,U)": "Gauss vanLeer;",  # "Gauss linear;",
        "div(phi,alpha)": "Gauss vanLeer01;",
        "div(((rho*nuEff)*dev2(T(grad(U)))))": "Gauss vanLeer;",
        "div(((rho*nuEff)*T(grad(U))))": "Gauss vanLeer;",
    }

    if use_turbulence:
        div_schemes["div(phi,U)"] = "Gauss limitedLinearV 1;"
        div_schemes["div(phi,k)"] = "Gauss limitedLinear 1;"
        div_schemes["div(phi,epsilon)"] = "Gauss limitedLinear 1;"
        div_schemes["div(phi,omega)"] = "Gauss limitedLinear 1;"
        div_schemes["div(phi,R)"] = "Gauss limitedLinear 1;"
        div_schemes["div(R)"] = "Gauss linear;"
        div_schemes["div(U)"] = "Gauss linear;"
        div_schemes["div(phi,nuTilda)"] = "Gauss limitedLinear 1;"

    s = ""
    for scheme in div_schemes:
        s += "\t" + scheme + "\t\t" + div_schemes[scheme] + "\n"

    with open(template_file, "r") as f:
        file_content = f.read()
    file_content = file_content.replace("<<<DIV_SCHEMES>>>", s)
    print("[init_psl] writing ...", path / "fvSchemes")
    with open(path / "fvSchemes", "w") as f:
        f.write(file_content)


def write_fv_solution(path, use_turbulence, template_file):
    solvers = {}

    if use_turbulence:
        solvers['"(k|epsilon|R|nuTilda|kFinal|epsilonFinal)"'] = {
            "solver": "smoothSolver;",
            "smoother": "GaussSeidel;",
            "tolerance": "1e-05;",
            "relTol": "0;",
        }

    s = ""
    for solver in solvers:
        s += "\t" + solver + "\n"
        s += "\t{\n"
        for parameter in solvers[solver]:
            s += "\t\t" + parameter + "\t\t" + solvers[solver][parameter] + "\n"
        s += "\t}\n"

    with open(template_file, "r") as f:
        file_content = f.read()
    file_content = file_content.replace("<<<SOLVERS>>>", s)
    print("[init_psl] writing ...", path / "fvSolution")
    with open(path / "fvSolution", "w") as f:
        f.write(file_content)


def write_decompose_par_dict(path, x, y, z, template_file):
    with open(template_file, "r") as f:
        file_content = f.read()
    file_content = file_content.replace(
        "<<<NUMBER_OF_SUBDOMAINS>>>", str(x * y * z)
    ).replace("<<<SUBDOMAINS>>>", "%d %d %d" % (x, y, z))

    print("[init_psl] writing ...", path / "decomposeParDict")
    with open(path / "decomposeParDict", "w") as f:
        f.write(file_content)


def write_psl_control(path, args):
    file_content = "// terrain patch\n"
    file_content += "pslTerrainPatchName terrain;\n"
    file_content += "// profile\n"
    file_content += "pslProfileFile stats;              // store mass, maxU per frame\n"
    if args.profile:
        file_content += "pslProfile 1;\n"
    else:
        file_content += "pslProfile 0;\n"
    file_content += "// entrainment method\n"
    file_content += (
        "pslEntrainmentMethod "
        + str(args.entrainment_method)
        + ";            // entrainment method [0 2]\n"
    )
    file_content += (
        "pslSnowCoverDensity 25;            // snow cover density [kg / m3]\n"
    )
    file_content += "pslUEntrainmentFactor " + str(args.u_factor) + ";\n"
    file_content += "pslAEntrainmentFactor " + str(args.a_factor) + ";\n"
    file_content += "pslFrontLengthFactor " + str(args.f_factor) + ";\n"
    print("[init_psl] writing ...", path / "pslControl")
    with open(path / "pslControl", "w") as f:
        f.write(file_content)


def write_all_run(path, set_fields, mesh_generator, template_file):
    with open(template_file, "r") as f:
        file_content = f.read()

    applications = "runApplication " + mesh_generator + "\n"
    if set_fields:
        applications += "runApplication setFields\n"

    file_content = file_content.replace("<<<APPLICATIONS>>>", str(applications))

    print("[init_psl] writing ...", path / os.path.basename(template_file))
    with open(path / os.path.basename(template_file), "w") as f:
        f.write(file_content)


def write_boundary_conditions(domain_type, use_turbulence, template_path, path, run_2d):
    conditions = {
        "alpha.snow": {},
        "U": {},
        "p_rgh": {},
    }

    if use_turbulence:
        conditions["k"] = {}
        conditions["epsilon"] = {}
        conditions["nut"] = {}
        conditions["nuTilda"] = {}

    for variable in conditions:
        conditions[variable]["top"] = {}
        conditions[variable]["terrain"] = {}
        conditions[variable]["walls"] = {}
        conditions[variable]["inlet"] = {}
        conditions[variable]["outlet"] = {}

    # input fields have only the wall option
    input_fields = [
        "alpha_inj",
        "alpha_inj_U",
        "dsl_U",
        "dsl_front_sdf",
        "dsl_h",
        "dsl_slope",
        "dsl_W",
        "dsl_divU",
        "Uinj",
        "gradAlpha1",
    ]
    for input_field in input_fields:
        conditions[input_field] = {}
        conditions[input_field]["walls"] = {}
        if run_2d:
            conditions[input_field]["walls"]["type"] = "empty"
        else:
            conditions[input_field]["walls"]["type"] = "zeroGradient"

    def createFreePatch(name):
        conditions["alpha.snow"][name]["type"] = "inletOutlet"
        conditions["alpha.snow"][name]["inletValue"] = "uniform 0"
        conditions["alpha.snow"][name]["value"] = "uniform 0"

        conditions["U"][name]["type"] = "pressureInletOutletVelocity"
        conditions["U"][name]["value"] = "uniform (0 0 0)"

        conditions["p_rgh"][name]["type"] = "totalPressure"
        conditions["p_rgh"][name]["p0"] = "uniform 0"

    def createWallPatch(name):
        conditions["alpha.snow"][name]["type"] = "fixedValue"
        conditions["alpha.snow"][name]["value"] = "uniform 0"

        conditions["U"][name]["type"] = "noSlip"

        conditions["p_rgh"][name]["type"] = "zeroGradient"

    def createOutletPatch(name):
        conditions["alpha.snow"][name]["type"] = "zeroGradient"

        conditions["U"][name]["type"] = "zeroGradient"

        conditions["p_rgh"][name]["type"] = "fixedValue"
        conditions["p_rgh"][name]["value"] = "uniform 0"

    def createInletPatch(name):
        conditions["alpha.snow"][name]["type"] = "fixedValue"
        conditions["alpha.snow"][name]["value"] = "uniform 0"

        conditions["U"][name]["type"] = "fixedValue"
        conditions["U"][name]["value"] = "uniform (-0.01 0 0)"

        conditions["p_rgh"][name]["type"] = "zeroGradient"

    if domain_type == "open":
        #
        #  . . . . . . . . . .
        #  .                 .
        #  .                 .
        #  ._________________.
        #

        createFreePatch("top")
        createFreePatch("inlet")
        createFreePatch("outlet")
        createFreePatch("walls")

        conditions["p_rgh"]["terrain"]["type"] = "fixedFluxPressure"
        conditions["p_rgh"]["terrain"]["value"] = "uniform 0"

    elif domain_type == "openBox":
        #
        #  . . . . . . . . . .
        #  |                 |
        #  |                 |
        #  |_________________|
        #
        createFreePatch("top")
        createWallPatch("inlet")
        createWallPatch("outlet")
        createWallPatch("walls")

        conditions["p_rgh"]["terrain"]["type"] = "fixedFluxPressure"
        conditions["p_rgh"]["terrain"]["value"] = "uniform 0"

    elif domain_type == "tunnel":
        #
        #   _________________
        #  .                 .
        #  .                 .
        #  ._________________.
        #
        createWallPatch("top")
        createFreePatch("inlet")
        createFreePatch("outlet")
        createWallPatch("walls")

        conditions["p_rgh"]["terrain"]["type"] = "fixedFluxPressure"
        conditions["p_rgh"]["terrain"]["value"] = "uniform 0"

    elif domain_type == "rev-wind-tunnel":
        #
        #   _________________
        #  <                 .
        #  <                 .
        #  <_________________.
        #
        createWallPatch("top")
        createInletPatch("inlet")
        createOutletPatch("outlet")
        createWallPatch("walls")

        conditions["p_rgh"]["terrain"]["type"] = "fixedFluxPressure"
        conditions["p_rgh"]["terrain"]["value"] = "uniform 0"

    elif domain_type == "rev-wind-open":
        #
        #  . . . . . . . . . .
        #  <                 .
        #  <                 .
        #  <_________________.
        #
        createFreePatch("top")
        createInletPatch("inlet")
        createOutletPatch("outlet")
        createWallPatch("walls")

        conditions["p_rgh"]["terrain"]["type"] = "fixedFluxPressure"
        conditions["p_rgh"]["terrain"]["value"] = "uniform 0"

    elif domain_type == "rev-full-wind-open":
        #
        #  . . . . . . . . . .
        #  <                 <
        #  <                 <
        #  <_________________<
        #
        createFreePatch("top")
        createInletPatch("inlet")
        createInletPatch("outlet")
        createWallPatch("walls")

        conditions["p_rgh"]["terrain"]["type"] = "fixedFluxPressure"
        conditions["p_rgh"]["terrain"]["value"] = "uniform 0"
    elif domain_type == "rev-full-wind-closed":
        #
        #   _________________
        #  <                 <
        #  <                 <
        #  <_________________<
        #
        createWallPatch("top")
        createInletPatch("inlet")
        createOutletPatch("outlet")
        createWallPatch("walls")

        conditions["p_rgh"]["terrain"]["type"] = "fixedFluxPressure"
        conditions["p_rgh"]["terrain"]["value"] = "uniform 0"

    elif domain_type == "L":
        #
        #  . . . . . . . . . .
        #  .                 |
        #  .                 |
        #  ._________________|
        #
        createWallPatch("top")
        createFreePatch("inlet")
        createWallPatch("outlet")
        createWallPatch("walls")

        conditions["p_rgh"]["terrain"]["type"] = "fixedFluxPressure"
        conditions["p_rgh"]["terrain"]["value"] = "uniform 0"

    elif domain_type == "wind-tunnel":
        conditions["alpha.snow"]["top"]["type"] = "inletOutlet"
        conditions["alpha.snow"]["top"]["inletValue"] = "uniform 0"
        conditions["alpha.snow"]["top"]["value"] = "uniform 0"
        conditions["alpha.snow"]["inlet"]["type"] = "inletOutlet"
        conditions["alpha.snow"]["inlet"]["inletValue"] = "uniform 0"
        conditions["alpha.snow"]["inlet"]["value"] = "uniform 0"
        conditions["alpha.snow"]["outlet"]["type"] = "inletOutlet"
        conditions["alpha.snow"]["outlet"]["inletValue"] = "uniform 0"
        conditions["alpha.snow"]["outlet"]["value"] = "uniform 0"
        conditions["alpha.snow"]["walls"]["type"] = "zeroGradient"

        conditions["U"]["walls"]["type"] = "slip"
        conditions["U"]["top"]["type"] = "pressureInletOutletVelocity"
        conditions["U"]["top"]["value"] = "uniform (0 0 0)"
        conditions["U"]["inlet"]["type"] = "pressureInletOutletVelocity"
        conditions["U"]["inlet"]["value"] = "uniform (0 0 0)"
        conditions["U"]["outlet"]["type"] = "pressureInletOutletVelocity"
        conditions["U"]["outlet"]["value"] = "uniform (0 0 0)"

        conditions["p_rgh"]["top"]["type"] = "totalPressure"
        conditions["p_rgh"]["top"]["p0"] = "uniform 0"
        conditions["p_rgh"]["inlet"]["type"] = "totalPressure"
        conditions["p_rgh"]["inlet"]["p0"] = "uniform 0"
        conditions["p_rgh"]["outlet"]["type"] = "totalPressure"
        conditions["p_rgh"]["outlet"]["p0"] = "uniform 0"
        conditions["p_rgh"]["walls"]["type"] = "fixedFluxPressure"
        conditions["p_rgh"]["walls"]["value"] = "uniform 0"
        conditions["p_rgh"]["terrain"]["type"] = "fixedFluxPressure"
        conditions["p_rgh"]["terrain"]["value"] = "uniform 0"
    elif domain_type == "bridge":
        conditions["alpha.snow"]["top"]["type"] = "inletOutlet"
        conditions["alpha.snow"]["top"]["inletValue"] = "uniform 0"
        conditions["alpha.snow"]["top"]["value"] = "uniform 0"
        conditions["alpha.snow"]["inlet"]["type"] = "inletOutlet"
        conditions["alpha.snow"]["inlet"]["inletValue"] = "uniform 0"
        conditions["alpha.snow"]["inlet"]["value"] = "uniform 0"
        conditions["alpha.snow"]["outlet"]["type"] = "inletOutlet"
        conditions["alpha.snow"]["outlet"]["inletValue"] = "uniform 0"
        conditions["alpha.snow"]["outlet"]["value"] = "uniform 0"
        conditions["alpha.snow"]["walls"]["type"] = "zeroGradient"

        conditions["U"]["walls"]["type"] = "slip"
        conditions["U"]["top"]["type"] = "noSlip"
        conditions["U"]["inlet"]["type"] = "pressureInletOutletVelocity"
        conditions["U"]["inlet"]["value"] = "uniform (0 0 0)"
        conditions["U"]["outlet"]["type"] = "pressureInletOutletVelocity"
        conditions["U"]["outlet"]["value"] = "uniform (0 0 0)"

        conditions["p_rgh"]["top"]["type"] = "zeroGradient"
        conditions["p_rgh"]["inlet"]["type"] = "totalPressure"
        conditions["p_rgh"]["inlet"]["p0"] = "uniform 0"
        conditions["p_rgh"]["outlet"]["type"] = "totalPressure"
        conditions["p_rgh"]["outlet"]["p0"] = "uniform 0"
        conditions["p_rgh"]["walls"]["type"] = "fixedFluxPressure"
        conditions["p_rgh"]["walls"]["value"] = "uniform 0"
        conditions["p_rgh"]["terrain"]["type"] = "zeroGradient"

        if use_turbulence:
            conditions["epsilon"]["top"]["type"] = "epsilonWallFunction"
            conditions["epsilon"]["top"]["value"] = "uniform 0.000765"
            conditions["epsilon"]["walls"]["type"] = "epsilonWallFunction"
            conditions["epsilon"]["walls"]["value"] = "uniform 0.000765"
            conditions["epsilon"]["terrain"]["type"] = "epsilonWallFunction"
            conditions["epsilon"]["terrain"]["value"] = "uniform 0.000765"
            conditions["epsilon"]["outlet"]["type"] = "zeroGradient"
            conditions["epsilon"]["inlet"]["type"] = "fixedValue"
            conditions["epsilon"]["inlet"]["value"] = "uniform 0.000765"

            conditions["k"]["top"]["type"] = "kqRWallFunction"
            conditions["k"]["top"]["value"] = "uniform 0.00325"
            conditions["k"]["walls"]["type"] = "kqRWallFunction"
            conditions["k"]["walls"]["value"] = "uniform 0.00325"
            conditions["k"]["terrain"]["type"] = "kqRWallFunction"
            conditions["k"]["terrain"]["value"] = "uniform 0.00325"
            conditions["k"]["outlet"]["type"] = "zeroGradient"
            conditions["k"]["inlet"]["type"] = "fixedValue"
            conditions["k"]["inlet"]["value"] = "uniform 0.00325"

            conditions["nut"]["top"]["type"] = "nutkWallFunction"
            conditions["nut"]["top"]["value"] = "uniform 0"
            conditions["nut"]["walls"]["type"] = "nutkWallFunction"
            conditions["nut"]["walls"]["value"] = "uniform 0"
            conditions["nut"]["terrain"]["type"] = "nutkWallFunction"
            conditions["nut"]["terrain"]["value"] = "uniform 0"
            conditions["nut"]["outlet"]["type"] = "zeroGradient"
            conditions["nut"]["inlet"]["type"] = "fixedValue"
            conditions["nut"]["inlet"]["value"] = "uniform 0"

            conditions["nuTilda"]["top"]["type"] = "zeroGradient"
            conditions["nuTilda"]["walls"]["type"] = "zeroGradient"
            conditions["nuTilda"]["terrain"]["type"] = "zeroGradient"
            conditions["nuTilda"]["outlet"]["type"] = "zeroGradient"
            conditions["nuTilda"]["inlet"]["type"] = "zeroGradient"

    if run_2d:
        for variable in conditions:
            conditions[variable]["walls"] = {}
            conditions[variable]["walls"]["type"] = "empty"

    for variable in conditions:
        template_file = template_path / variable
        with open(template_file, "r") as f:
            file_content = f.read()

        for patch in conditions[variable]:
            condition = ""
            for entry in conditions[variable][patch]:
                condition += (
                    "\t\t\t"
                    + entry
                    + "\t\t\t"
                    + conditions[variable][patch][entry]
                    + ";\n"
                )
            file_content = file_content.replace(
                "<<<" + patch.upper() + "_CONDITIONS>>>", condition
            )

        with open(path / variable, "w") as f:
            f.write(file_content)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", type=Path, help="sim path")
    parser.add_argument("--template-dir", type=Path, help="template files dir")
    parser.add_argument(
        "--write-format",
        type=str,
        help="output file type (ascii|binary)",
        default="ascii",
    )
    parser.add_argument("--dt", type=float, help="time step", default=0.0001)
    parser.add_argument("--wt", type=float, help="write time", default=0.05)
    parser.add_argument("--et", type=float, help="end time", default=10)
    parser.add_argument("--px", type=float, help="x domains", default=5)
    parser.add_argument("--py", type=float, help="y domains", default=4)
    parser.add_argument("--pz", type=float, help="z domains", default=1)
    parser.add_argument(
        "--snow-density", type=float, help="snow cloud density kg/m3", default=7
    )
    parser.add_argument(
        "--air-density", type=float, help="air density kg/m3", default=1.225
    )
    parser.add_argument("--set-fields", action="store_true")
    parser.add_argument(
        "--mesh-generator",
        type=str,
        help="blockMesh|pMesh|cartesianMesh|makeFaMesh|slopeMesh",
        default="blockMesh",
    )
    parser.add_argument("--max-co", type=str, default="0.2")
    parser.add_argument("--dab", type=str, default="2e-04")
    parser.add_argument("--snow-viscosity", type=str, default="1e-04")
    parser.add_argument(
        "--domain-type", type=str, help="openBox | tunnel", default="openBox"
    )
    parser.add_argument("--run-2d", action="store_true")
    parser.add_argument("--use-turbulence", action="store_true")
    parser.add_argument("--entrainment-method", type=int, help="[0 2]", default=2)
    parser.add_argument(
        "--u-factor", type=float, help="velocity entrainment factor", default=1
    )
    parser.add_argument(
        "--a-factor", type=float, help="mass entrainment factor", default=1
    )
    parser.add_argument(
        "--f-factor", type=float, help="front length factor", default=10
    )
    parser.add_argument("--profile", action="store_true")
    args = parser.parse_args()

    if args.o.exists() and os.listdir(args.o):
        print("[init_psl] output folder not empty")
        print("[init_psl][output folder] " + str(args.o))
        exit(1)
    elif not args.o.exists():
        os.makedirs(args.o)

    constant_dir = args.o / "constant"
    system_dir = args.o / "system"
    variables_dir = args.o / "0.orig"

    # copy template files into new directory
    shutil.copytree(
        args.template_dir,
        args.o,
        symlinks=False,
        ignore=None,
        ignore_dangling_symlinks=False,
        dirs_exist_ok=True,
    )
    if not args.use_turbulence:
        os.remove(variables_dir / "k")
        os.remove(variables_dir / "epsilon")
        os.remove(variables_dir / "nut")
        os.remove(variables_dir / "nuTilda")

    write_psl_control(system_dir, args)
    write_control_dict_file(
        system_dir,
        args.dt,
        args.wt,
        args.et,
        args.write_format,
        args.max_co,
        args.template_dir / "system/controlDict",
    )
    write_transport_properties(
        constant_dir,
        args.snow_density,
        args.air_density,
        args.template_dir / "constant/transportProperties",
        args.dab,
        args.snow_viscosity,
    )
    write_turbulence_properties(
        constant_dir,
        args.use_turbulence,
        args.template_dir / "constant/turbulenceProperties",
    )
    write_fv_schemes(
        system_dir, args.use_turbulence, args.template_dir / "system/fvSchemes"
    )
    write_fv_solution(
        system_dir, args.use_turbulence, args.template_dir / "system/fvSolution"
    )
    write_decompose_par_dict(
        system_dir,
        args.px,
        args.py,
        args.pz,
        args.template_dir / "system/decomposeParDict",
    )
    write_all_run(
        args.o, args.set_fields, args.mesh_generator, args.template_dir / "Allrun"
    )
    write_all_run(
        args.o,
        args.set_fields,
        args.mesh_generator,
        args.template_dir / "Allrun-parallel",
    )
    write_boundary_conditions(
        args.domain_type,
        args.use_turbulence,
        args.template_dir / "0.orig",
        variables_dir,
        args.run_2d,
    )
