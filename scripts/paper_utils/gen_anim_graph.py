import os
import argparse
import shutil
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.animation as animation


def plumes(frame, roots, args):
    df_plumes = pd.DataFrame()

    for root in roots:
        name = root.split("/")[-args.no]

        plumes_filename = root + "/stats/" + str(frame).zfill(6) + "_plumes.csv"
        labels = ["distance", "height"]
        data = pd.read_csv(
            plumes_filename, sep=" ", names=labels, header=None, index_col=False
        )
        data["name"] = name
        df_plumes = pd.concat([df_plumes, data], axis=0)

    return df_plumes


fig, ax = plt.subplots()
roots = []

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=Path, help="stats directory")
    parser.add_argument("-fs", type=int, help="frame count")
    parser.add_argument("--d2", action="store_true")
    parser.add_argument("--no", type=int, default=3)
    args = parser.parse_args()

    # go over all directories and find stats folders
    for dirpath, dirs, files in os.walk(args.i):
        # check for stats dir
        if "stats" in dirs:
            roots.append(dirpath)

    plt.title("Plume Profile", y=0.95)
    plt.xlim(0, 800)
    plt.xlim(100, 800)

    def update(frame):
        if frame == 0:
            return
        global ax
        # ax.cla()
        df_plumes = plumes(frame, roots, args)

        ax = sns.relplot(
            data=df_plumes, kind="line", x="distance", y="height", hue="name"
        )
        ax.set(ylabel="Height [m]", xlabel="Position [m]", title="")
        sns.move_legend(
            ax, "upper center", bbox_to_anchor=(0.22, 0.95), title=None, frameon=True
        )

    ani = animation.FuncAnimation(fig, update, frames=args.fs)
    ani.save("test.mp4", dpi=200, fps=20)
