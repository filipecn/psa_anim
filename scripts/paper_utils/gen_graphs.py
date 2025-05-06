import os
import argparse
import shutil
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re


def processLog(logfile):
    with open(logfile, "r") as f:
        lines = f.readlines()
        deltaTCount = 1
        courantCount = 1
        executionCount = 1
        deltaTacc = 0
        courantAcc = 0
        courantMaxAcc = 0
        lastExecution = 0
        prefix_delta = "deltaT = "
        prefix_courant = "Courant Number mean:"
        prefix_exec = "ExecutionTime = "
        courant_number_means = []
        max_courant_number_means = []
        deltaTs = []
        reg = r"[-+]?(?:\d*\.*\d*[e]?[-+]?\d+)"
        for line in lines:
            if prefix_delta in line:
                r = re.findall(reg, line)
                deltaTacc += float(r[0])
                deltaTCount += 1
                deltaTs.append(float(r[0]))
            elif prefix_courant in line:
                r = re.findall(reg, line)
                courantAcc += float(r[0])
                courantMaxAcc += float(r[1])
                courantCount += 1
                courant_number_means.append(float(r[0]))
                max_courant_number_means.append(float(r[1]))
            elif prefix_exec in line:
                r = re.findall(reg, line)
                lastExecution = float(r[0])
                executionCount += 1
        print(logfile)
        print(
            deltaTacc / deltaTCount,
            courantAcc / courantCount,
            courantMaxAcc / courantCount,
            lastExecution,
        )
        return [courant_number_means, max_courant_number_means, deltaTs, lastExecution]


def merge_stats(folder, out):
    nthreads = 20
    data = []
    for thread in range(nthreads):
        with open(folder + "/" + "stats" + str(thread)) as file:
            lines = file.readlines()
            i = 0
            for line in lines:
                values = [float(x) for x in line.split()]
                if len(data) == i:
                    data.append(values)
                else:
                    # 0 time
                    data[i][0] == values[0]
                    # 1  total_mass
                    data[i][1] += values[1]
                    # 2  total_entrained_mass
                    data[i][2] += values[2]
                    # 3  min_mag_u
                    data[i][3] = min(data[i][3], values[3])
                    # 4  max_mag_u
                    data[i][4] = max(data[i][4], values[4])
                    # 5  total_kinetic_energy
                    data[i][5] += values[5]
                    # 6  total_entrained_kinetic_energy
                    data[i][6] += values[6]
                    # 7  total_potential_energy
                    data[i][7] += values[7]
                    # 8  max_alpha_inj
                    data[i][8] = max(data[i][8], values[8])
                    # 9  max_mag_alpha_inj_U
                    data[i][9] = max(data[i][9], values[9])
                    # 10 dsl_total_mass
                    data[i][10] += values[10]
                i += 1
    with open(out, "w") as file:
        for d in data:
            file.write(" ".join([str(x) for x in d]) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=Path, help="stats directory")
    parser.add_argument("-f", type=int, help="frame number")
    parser.add_argument("--d2", action="store_true")
    parser.add_argument("--no", type=int, default=3)
    parser.add_argument("--names", type=str, default="")
    parser.add_argument("--max-frames", type=int, default=500)
    args = parser.parse_args()

    roots = []
    # go over all directories and find stats folders
    for dirpath, dirs, files in os.walk(args.i):
        # check for stats dir
        if "stats" in dirs:
            roots.append(dirpath)

    df_plumes = pd.DataFrame()
    df_front = pd.DataFrame()
    df_stats = pd.DataFrame()

    dvalues = ["0", "1", "2", "3"]

    if args.names == "densities":
        dvalues = ["1.4", "2.5", "5.0", "7.0"]

    if args.names == "a":
        dvalues = ["0.01", "0.1", "0.5", "0.9"]

    if args.names == "u":
        dvalues = ["1.0", "2.0", "4.0", "8.0"]

    if args.names == "e":
        dvalues = ["15", "50", "100", "400"]

    if args.names == "lfront":
        dvalues = ["10", "20", "40", "80"]

    if args.names == "convergence":
        dvalues = ["0.1", "0.5", "1.0", "5.0", "10.0"]

    dnames = {}
    i = 0
    for v in dvalues:
        dnames["psl_" + str(i)] = dvalues[i]
        i += 1

    i = 0
    delta_t_data = []
    max_courant_number_data = []
    courant_number_data = []
    duration_data = []
    for root in roots:
        i += 1
        name = root.split("/")[-args.no]
        if name in dnames:
            name = dnames[name]

        # stats
        psl_folder = "psl"
        if args.d2:
            psl_folder = "psl2"

        stats_labels_filename = "/".join([root, psl_folder, "stats_labels"])

        log_filename = "/".join([root, psl_folder, "log.pslFoam"])
        log_data = processLog(log_filename)
        delta_t_data.append(log_data[2])
        max_courant_number_data.append(log_data[1])
        courant_number_data.append(log_data[0])
        duration_data.append(log_data[3])

        # if there is more than one stats, combine all
        stats_filename = "/".join([root, psl_folder, "stats0"])
        if os.path.exists("/".join([root, psl_folder, "stats1"])):
            stats_filename = "merged_stats" + str(i)
            merge_stats(root + "/" + psl_folder, stats_filename)

        labels = []
        with open(stats_labels_filename, "r") as f:
            labels = f.readline().split()
        data = pd.read_csv(
            stats_filename, sep=" ", names=labels, header=None, index_col=False
        )
        data["name"] = name
        df_stats = pd.concat([df_stats, data.head(args.max_frames)], axis=0)

        # plumes
        plumes_filename = root + "/stats/" + str(args.f).zfill(6) + "_plumes.csv"
        labels = ["distance", "height"]
        data = pd.read_csv(
            plumes_filename, sep=" ", names=labels, header=None, index_col=False
        )
        data["name"] = name
        df_plumes = pd.concat([df_plumes, data], axis=0)

        # front
        front_filename = root + "/stats/front_data.csv"
        labels = ["time", "front position", "front velocity"]
        data = pd.read_csv(
            front_filename, sep=",", names=labels, header=None, index_col=False
        )
        # data["front position"] = [
        #    x - data["front position"][0] for x in data["front position"]
        # ]
        data["name"] = name
        df_front = pd.concat([df_front, data.head(args.max_frames)], axis=0)

    df_stats["total_volume"] = df_stats["total_mass"] / 1.4

    sns.set_context(
        "paper", rc={"font.size": 12, "axes.titlesize": 12, "axes.labelsize": 12}
    )
    # stats
    ax = sns.relplot(data=df_stats, kind="line", x="time", y="total_mass", hue="name")
    ax.set(ylabel="Total Mass [kg]", xlabel="Time [s]", title="")
    plt.title("Total Mass", y=0.95)
    sns.move_legend(
        ax, "upper center", bbox_to_anchor=(0.25, 0.95), title=None, frameon=True
    )
    plt.savefig("total_mass.png")

    ax = sns.relplot(data=df_stats, kind="line", x="time", y="total_volume", hue="name")
    ax.set(ylabel="Total Volume [m3]", xlabel="Time [s]", title="")
    plt.title("Total Volume", y=0.95)
    sns.move_legend(
        ax, "upper center", bbox_to_anchor=(0.25, 0.95), title=None, frameon=True
    )
    plt.savefig("total_volume.pdf")

    ax = sns.relplot(
        data=df_stats, kind="line", x="time", y="dsl_total_mass", hue="name"
    )
    ax.set(ylabel="Total Mass [kg]", xlabel="Time [s]", title="")
    plt.title("Total Mass", y=0.95)
    sns.move_legend(
        ax, "upper center", bbox_to_anchor=(0.25, 0.95), title=None, frameon=True
    )
    plt.savefig("dsl_total_mass.png")

    # max u
    ax = sns.relplot(data=df_stats, kind="line", x="time", y="max_mag_u", hue="name")
    ax.set(ylabel="Velocity [m/s]", xlabel="Time [s]", title="")
    plt.title("Maximum Velocity", y=0.95)
    sns.move_legend(
        ax, "upper center", bbox_to_anchor=(0.25, 0.95), title=None, frameon=True
    )
    plt.savefig("max_mag_u.png")

    # front
    ax = sns.relplot(
        data=df_front, kind="line", x="time", y="front position", hue="name"
    )
    ax.set(ylabel="Distance [m]", xlabel="time [s]", title="")
    plt.title("Front Position", y=0.95)
    sns.move_legend(
        ax, "upper center", bbox_to_anchor=(0.25, 0.95), title=None, frameon=True
    )
    plt.savefig("frontpos.pdf", format="pdf")

    ax = sns.relplot(
        data=df_front, kind="line", x="time", y="front velocity", hue="name"
    )
    ax.set(ylabel="Velocity [m/s]", xlabel="time [s]", title="")
    sns.move_legend(
        ax, "upper center", bbox_to_anchor=(0.22, 0.95), title=None, frameon=True
    )
    plt.title("Front Velocity", y=0.95)
    plt.savefig("frontvel.pdf", format="pdf")

    # plumes
    ax = sns.relplot(data=df_plumes, kind="line", x="distance", y="height", hue="name")
    ax.set(ylabel="Height [m]", xlabel="Position [m]", title="")
    sns.move_legend(
        ax, "upper center", bbox_to_anchor=(0.22, 0.95), title=None, frameon=True
    )
    plt.title("Plume Profile", y=0.95)
    plt.xlim(0, 800)
    plt.xlim(0, 800)
    plt.savefig("plumes.pdf", format="pdf")

    # ax = sns.boxplot(y=courant_number_data)
    # plt.savefig("courant.png")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.title("Mean Courant Number")
    ax.set_xticklabels(dvalues)
    plt.boxplot(courant_number_data)
    plt.savefig("courant.pdf", format="pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.title("Time Step")
    ax.set_xticklabels(dvalues)
    plt.boxplot(delta_t_data)
    plt.savefig("dt.pdf", format="pdf")

    print(duration_data)
    print(sum(duration_data) / len(duration_data))
