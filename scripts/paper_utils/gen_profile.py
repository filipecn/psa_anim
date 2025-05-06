import os
import argparse
import shutil
from pathlib import Path
import re


def readPSLCellCount(path: Path):
    if not path.exists():
        return 0
    with open(path, "r") as file:
        for line in file:
            matches = re.findall(r"nCells: (\d+)", line)
            if len(matches):
                return int(matches[-1])
            matches = re.findall(r"cells: *(\d+)", line)
            if len(matches):
                return int(matches[-1])
    return 0


def readDSLCellCount(path: Path):
    if not path.exists():
        return 0
    with open(path, "r") as file:
        for line in file:
            matches = re.findall(r"nFaces: (\d+)", line)
            if len(matches):
                return int(matches[-1])
            matches = re.findall(r"Number of faces: (\d+)", line)
            if len(matches):
                return int(matches[-1])
    return 0


def readTotalTime(path: Path):
    if not path.exists():
        return 0
    found = None
    with open(path, "r") as file:
        for line in file:
            matches = re.findall(r"ExecutionTime = (\d*.\d*)", line)
            if len(matches):
                found = matches
    if found:
        return float(found[0])
    return 0


def readMainLoopIterations(path: Path):
    if not path.exists():
        return 0
    found = 0
    with open(path, "r") as file:
        for line in file:
            if "ExecutionTime =" in line:
                found += 1
    return found


def readProfileTimes(path: Path):
    if not path.exists():
        return {}
    times = {}
    with open(path, "r") as file:
        for line in file:
            txt = line.split()
            name = txt[0]
            time = float(txt[1])
            count = int(txt[2])
            if name not in times:
                times[name] = {"time": 0, "count": 0}
            times[name]["time"] += time
            times[name]["count"] += count

    # for item in times:
    #    prof_item = times[item]
    #    print(
    #        "{:30s} {:5.4f}% {:10.4f} {}".format(
    #            item,
    #            prof_item["time"] / times["totalTime"]["time"] * 100,
    #            prof_item["time"],
    #            prof_item["count"],
    #        )
    #    )
    return times


def readPSLParameters(path: Path):
    parameters = {}
    if not path.exists():
        return parameters
    with open(path, "r") as file:
        for line in file:
            matches = re.findall(r"([a-zA-Z]+) ([0-9\.]+);", line)
            if len(matches):
                parameters[matches[0][0]] = str(matches[0][1])
    return parameters


def readDSLTransportProperties(path: Path):
    parameters = {}
    if not path.exists():
        return parameters
    with open(path, "r") as file:
        for line in file:
            if line[:2] == "//":
                continue
            matches = re.findall(
                r"[ ]*([a-zA-Z0-9]+)[ ]+([a-zA-Z0-9]+)[ ]+\[([ \-0-9]+)\][ ]+([\-e0-9\.]+);",
                line,
            )
            if len(matches):
                parameters[matches[0][0]] = str(matches[0][-1])
    return parameters


def readPSLTransportProperties(path: Path):
    parameters = {}
    if not path.exists():
        return parameters
    with open(path, "r") as file:
        for line in file:
            if line[:2] == "//":
                continue
            matches = re.findall(
                r"[ ]*([a-zA-Z0-9]+)[ ]+([\-e0-9\.]+);",
                line,
            )
            if len(matches):
                name = matches[0][0]
                if name in parameters:
                    name += "_air"
                parameters[name] = str(matches[0][-1])
    return parameters


def dslMesh(path: Path):
    sizes = {}
    if not path.exists():
        return sizes
    with open(path, "r") as file:
        for line in file:
            matches = re.findall(
                r"Number of faces: ([0-9]+)",
                line,
            )
            if len(matches):
                sizes["faces"] = str(matches[0])
                continue
            matches = re.findall(
                r"[ ]*min = ([0-9\.e]+) max = ([0-9\.e]+)",
                line,
            )
            if len(matches) and "min_area" not in sizes:
                sizes["min_area"] = str(matches[0][0])
                sizes["max_area"] = str(matches[0][1])
    return sizes


def pslMesh(path: Path):
    sizes = {}
    if not path.exists():
        return sizes
    with open(path, "r") as file:
        for line in file:
            matches = re.findall(
                r"[' ']+cells:[ ]+([0-9]+)",
                line,
            )
            if len(matches):
                sizes["cells"] = str(matches[0])
            matches = re.findall(
                r"([a-zA-Z]+) volume = ([0-9\.e+]+)",
                line,
            )
            for match in matches:
                if match[1][-1] == ".":
                    sizes[match[0] + " volume"] = match[1][:-1]
                else:
                    sizes[match[0] + " volume"] = match[1]
    sizes["avg"] = str(
        (float(sizes["Total volume"]) / float(sizes["cells"])) ** (1.0 / 3.0)
    )
    return sizes


def readMaxCo(path: Path):
    if not path.exists():
        return 0.2

    with open(path, "r") as file:
        for line in file:
            if "maxCo" in line:
                return float(line.strip().strip(";").split()[1])
    return 2.0


def readSolutionParameters(path: Path):
    parms = {"npiso": 2, "ncorr": 2}
    if not path.exists():
        return parms

    with open(path, "r") as file:
        for line in file:
            if "nCorrectors" in line:
                parms["npiso"] = int(line.strip().strip(";").split()[1])
            if "nNonOrthogonalCorrectors" in line:
                parms["ncorr"] = int(line.strip().strip(";").split()[1])
    return parms


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=Path, help="root directory")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--table", action="store_true")
    parser.add_argument("--profile", action="store_true")
    args = parser.parse_args()

    profile = {}

    for sim_name in os.listdir(args.i):
        sim_path = args.i / sim_name
        profile[sim_name] = {}
        profile[sim_name]["profile"] = readProfileTimes(sim_path / "profile_times_0")

        # sum run times for psl's
        profile[sim_name]["psl_total_time"] = 0
        profile[sim_name]["main loop iterations"] = 0
        for file in [x for x in os.listdir(sim_path) if "log.pslFoam" in x]:
            profile[sim_name]["psl_total_time"] += readTotalTime(sim_path / file)
            profile[sim_name]["main loop iterations"] += readMainLoopIterations(
                sim_path / file
            )

        profile[sim_name]["psl_cell_count"] = readPSLCellCount(sim_path / "checkMesh")
        profile[sim_name]["dsl_cell_count"] = readDSLCellCount(sim_path / "checkFaMesh")
        profile[sim_name]["dsl_total_time"] = readTotalTime(
            sim_path / "log.faSavageHutterFoam"
        )
        profile[sim_name]["parameters"] = readPSLParameters(sim_path / "pslControl")
        profile[sim_name]["dsl_transport"] = readDSLTransportProperties(
            sim_path / "dslTransportProperties"
        )
        profile[sim_name]["psl_transport"] = readPSLTransportProperties(
            sim_path / "pslTransportProperties"
        )
        profile[sim_name]["dsl_mesh"] = dslMesh(sim_path / "checkFaMesh")
        profile[sim_name]["psl_mesh"] = pslMesh(sim_path / "checkMesh")
        profile[sim_name]["max co"] = readMaxCo(sim_path / "controlDict")
        profile[sim_name]["solution"] = readSolutionParameters(sim_path / "fvSolution")

    # print parameters
    if args.all:
        parameter_names = []
        for sim in profile:
            prof = profile[sim]
            for item in prof["parameters"]:
                if item != "pslProfile":
                    parameter_names.append(item)
            break
        s = "{:20s}".format("name")
        s += " ".join(["{:15s}".format(x) for x in parameter_names])
        print(s)
        for sim in profile:
            s = "{:20s}".format(sim)
            prof = profile[sim]["parameters"]
            s += " ".join(["{:20s}".format(prof[x]) for x in parameter_names])
            print(s)

        # print transport
        parameter_names = []
        for sim in profile:
            prof = profile[sim]
            for item in prof["dsl_transport"]:
                parameter_names.append(item)
            break
        s = "{:15s}".format("name")
        s += " ".join(["{:15s}".format(x) for x in parameter_names])
        print(s)
        for sim in profile:
            s = "{:15s}".format(sim)
            prof = profile[sim]["dsl_transport"]
            s += " ".join(["{:15s}".format(prof[x]) for x in parameter_names])
            print(s)

        # print transport
        parameter_names = []
        for sim in profile:
            prof = profile[sim]
            for item in prof["psl_transport"]:
                if item != "version":
                    parameter_names.append(item)
            break
        s = "{:15s}".format("name")
        s += " ".join(["{:15s}".format(x) for x in parameter_names])
        print(s)
        for sim in profile:
            s = "{:15s}".format(sim)
            prof = profile[sim]["psl_transport"]
            s += " ".join(["{:15s}".format(prof[x]) for x in parameter_names])
            print(s)

        # dsl mesh
        s = "{:15s}".format("name")
        parameter_names = []
        for sim in profile:
            prof = profile[sim]
            for item in prof["dsl_mesh"]:
                parameter_names.append(item)
            break
        s = "{:15s}".format("name")
        s += " ".join(["{:15s}".format(x) for x in parameter_names])
        print(s)
        for sim in profile:
            s = "{:15s}".format(sim)
            prof = profile[sim]["dsl_mesh"]
            s += " ".join(["{:15s}".format(prof[x]) for x in parameter_names])
            print(s)

        # psl mesh
        s = "{:15s}".format("name")
        parameter_names = []
        for sim in profile:
            prof = profile[sim]
            for item in prof["psl_mesh"]:
                parameter_names.append(item)
            break
        s = "{:15s}".format("name")
        s += " ".join(["{:15s}".format(x) for x in parameter_names])
        print(s)
        for sim in profile:
            s = "{:15s}".format(sim)
            prof = profile[sim]["psl_mesh"]
            s += " ".join(["{:15s}".format(prof[x]) for x in parameter_names])
            print(s)

        print("main loop iterations from log:")
        for sim in profile:
            s = "{:15s}".format(sim)
            s += str(profile[sim]["main loop iterations"])
            print(s)

    profile["alpine"]["duration"] = 42
    profile["wolfsgrube"]["duration"] = 78
    profile["logo"]["duration"] = 40
    profile["slope25"]["duration"] = 30
    profile["slope30"]["duration"] = 30
    profile["slope35"]["duration"] = 30
    profile["ramp"]["duration"] = 35
    profile["river"]["duration"] = 35

    profile["niobe_open"]["duration"] = 100
    profile["niobe_closed"]["duration"] = 100
    profile["logo_original"]["duration"] = 90
    profile["logo_optimized"]["duration"] = 90
    profile["slope_low"]["duration"] = 30
    profile["slope_high"]["duration"] = 30

    for sim in profile:
        if "slope" in sim:
            profile[sim]["dim"] = 2
        else:
            profile[sim]["dim"] = 3

    # print profile

    if args.profile or args.all:

        def print_profile(name):
            print(name)
            print("=" * (40 + 10 + 20 + 10 + 3 * 3 + 1))
            print(
                "{:40s} | {:10s}  | {:20s} | {:10s}".format(
                    "section", "%", "total time", "calls"
                )
            )
            prof = profile[name]
            for item in prof["profile"]:
                prof_item = prof["profile"][item]
                print(
                    "{:40s} | {:10.4f}% | {:20.4f} | {:10d}".format(
                        item,
                        prof_item["time"] / prof["profile"]["totalTime"]["time"] * 100,
                        prof_item["time"],
                        prof_item["count"],
                    )
                )
            print("\n\n")

        print_profile("logo_original")
        print_profile("logo_optimized")
        print_profile("slope_high")
        print_profile("slope_low")
        print_profile("niobe_open")
        print_profile("niobe_closed")

    if args.table or args.all:
        table_names = {
            "scene": "{:30s} ",
            "dim": "& {:3s} ",
            "dur.": "& {:4s} ",
            "# DSL": "& {:10s} ",
            "# PSL": "& {:10s} ",
            "dx DSL": "& {:10s} ",
            "dx PSL": "& {:10s} ",
            "npiso": "& {:5s} ",
            "ncorr": "& {:5s} ",
            "co max": "& {:6s} ",
            "per dt": "& {:10s} ",
            "t DSL": "& {:10s} ",
            "t PSL": "& {:10s} ",
            "t tot": "& {:10s} ",
        }

        table_columns = {
            "scene": "{:30s} ",
            "dim": "& {:3d} ",
            "dur.": "& {:4d} ",
            "# DSL": "& {:10d} ",
            "# PSL": "& {:10d} ",
            "dx DSL": "& {:10.2f} ",
            "dx PSL": "& {:10.2f} ",
            "npiso": "& {:5d} ",
            "ncorr": "& {:5d} ",
            "co max": "& {:6.2f} ",
            "per dt": "& {:10.2f} ",
            "t DSL": "& {:10.2f} ",
            "t PSL": "& {:10.2f} ",
            "t tot": "& {:10.2f} ",
        }

        table_columns_names_txt = ""
        table_columns_names_txt += table_names["scene"]
        table_columns_names_txt += table_names["dim"]
        table_columns_names_txt += table_names["dur."]
        table_columns_names_txt += table_names["# DSL"]
        table_columns_names_txt += table_names["# PSL"]
        table_columns_names_txt += table_names["dx DSL"]
        table_columns_names_txt += table_names["dx PSL"]
        table_columns_names_txt += table_names["npiso"]
        table_columns_names_txt += table_names["ncorr"]
        table_columns_names_txt += table_names["co max"]
        table_columns_names_txt += table_names["per dt"]
        table_columns_names_txt += table_names["t DSL"]
        table_columns_names_txt += table_names["t PSL"]
        table_columns_names_txt += table_names["t tot"]

        print(
            table_columns_names_txt.format(
                "scene",
                "dim",
                "dur.",
                "# DSL",
                "# PSL",
                "dx DSL",
                "dx PSL",
                "npiso",
                "ncorr",
                "co max",
                "per dt",
                "t DSL",
                "t PSL",
                "t tot",
            )
        )

        table_columns_txt = ""
        table_columns_txt += table_columns["scene"]
        table_columns_txt += table_columns["dim"]
        table_columns_txt += table_columns["dur."]
        table_columns_txt += table_columns["# DSL"]
        table_columns_txt += table_columns["# PSL"]
        table_columns_txt += table_columns["dx DSL"]
        table_columns_txt += table_columns["dx PSL"]
        table_columns_txt += table_columns["npiso"]
        table_columns_txt += table_columns["ncorr"]
        table_columns_txt += table_columns["co max"]
        table_columns_txt += table_columns["per dt"]
        table_columns_txt += table_columns["t DSL"]
        table_columns_txt += table_columns["t PSL"]
        table_columns_txt += table_columns["t tot"]

        def print_paper_table_line(name):
            prof = profile[name]
            # compute per dt time
            per_dt = -1.0
            if "mainLoopIterationTime" in prof["profile"]:
                per_dt = (
                    prof["profile"]["mainLoopIterationTime"]["time"]
                    / prof["profile"]["mainLoopIterationTime"]["count"]
                )
            else:
                # writting output files takes approximately 4%
                per_dt = 0.96 * (prof["psl_total_time"] / prof["main loop iterations"])
            print(
                table_columns_txt.format(
                    name,
                    prof["dim"],
                    prof["duration"],
                    prof["dsl_cell_count"],
                    prof["psl_cell_count"],
                    float(prof["dsl_mesh"]["max_area"]),
                    float(prof["psl_mesh"]["avg"]),
                    prof["solution"]["npiso"],
                    prof["solution"]["ncorr"],
                    prof["max co"],
                    per_dt,
                    prof["dsl_total_time"],
                    prof["psl_total_time"],
                    prof["dsl_total_time"] + prof["psl_total_time"],
                )
            )

        print_paper_table_line("alpine")
        print_paper_table_line("wolfsgrube")
        print_paper_table_line("river")
        print_paper_table_line("ramp")
        print_paper_table_line("slope25")
        print_paper_table_line("slope30")
        print_paper_table_line("slope35")
        print_paper_table_line("logo_original")
        print_paper_table_line("logo_optimized")
        print_paper_table_line("niobe_closed")
        print_paper_table_line("niobe_open")
        print_paper_table_line("slope_low")
        print_paper_table_line("slope_high")
