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
    return 0


def readDSLCellCount(path: Path):
    if not path.exists():
        return 0
    with open(path, "r") as file:
        for line in file:
            matches = re.findall(r"nFaces: (\d+)", line)
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

    print(times)
    for item in times:
        prof_item = times[item]
        print(
            "{:30s} {:5.4f}% {:10.4f} {}".format(
                item,
                prof_item["time"] / times["totalTime"]["time"] * 100,
                prof_item["time"],
                prof_item["count"],
            )
        )
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=Path, help="root director")
    args = parser.parse_args()

    profile = {}

    for sim_name in os.listdir(args.i):
        sim_path = args.i / sim_name
        profile[sim_name] = {}
        profile[sim_name]["profile"] = readProfileTimes(sim_path / "profile_times_0")

        # sum run times for psl's
        profile[sim_name]["psl_total_time"] = 0
        for file in [x for x in os.listdir(sim_path) if "log.pslFoam" in x]:
            profile[sim_name]["psl_total_time"] += readTotalTime(sim_path / file)

        print(profile[sim_name]["psl_total_time"])

        profile[sim_name]["psl_cell_count"] = readPSLCellCount(
            sim_path / "log.blockMesh"
        )
        profile[sim_name]["dsl_cell_count"] = readDSLCellCount(
            sim_path / "log.makeFaMesh"
        )
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

    profile["alpine"]["duration"] = 42
    profile["wolfsgrube"]["duration"] = 78
    profile["logo"]["duration"] = 40
    profile["slope25"]["duration"] = 30
    profile["slope30"]["duration"] = 30
    profile["slope35"]["duration"] = 30
    profile["ramp"]["duration"] = 35
    profile["river"]["duration"] = 35

    # print parameters
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

    # print profile
    print(
        "{:10s} & {:10s} & {:10s} & {:10s} & {:10s} & {:10s} & {:10s}".format(
            "scene", "duration", "DSL", "PSL", "total", "DSL", "PSL"
        )
    )
    for sim in profile:
        prof = profile[sim]
        print(
            "{:10s} & {:10.2f} & {:10d} & {:10d} & {:10.2f} & {:10.2f} & {:10.2f}".format(
                sim,
                prof["duration"],
                prof["dsl_cell_count"],
                prof["psl_cell_count"],
                prof["dsl_total_time"] + prof["psl_total_time"],
                prof["dsl_total_time"],
                prof["psl_total_time"],
            )
        )

    # print profile
    prof = profile["logo"]
    for item in prof["profile"]:
        prof_item = prof["profile"][item]
        print(
            "{:30s} {:5.4f}% {:10.4f} {}".format(
                item,
                prof_item["time"] / prof["profile"]["totalTime"]["time"] * 100,
                prof_item["time"],
                prof_item["count"],
            )
        )
