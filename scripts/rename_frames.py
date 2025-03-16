# this script rename simulation frame files
# openfoam name frames with time (a float), this script simply multiply the float
# so we can easily sort the files later

import argparse
from pathlib import Path
import os
import shutil

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', type=Path, help="input dir")
    parser.add_argument('--o', type=Path, help="output dir")
    args = parser.parse_args()

    files = [x for x in os.listdir(args.i) if (args.i / x).is_file() and x[-3:] == "png"]
    for file in files:
        try:
            f = float(file[:-4])
        except ValueError:
            continue
        n = int(f * 100)
        name = f"{n:06d}" + file[-4:]
        shutil.copy(args.i / file, args.o / name)
