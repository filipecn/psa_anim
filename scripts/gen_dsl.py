#!/bin/python3
# this script generates a procedural dsl output data set of some types of motion:
# - rotating disk
# on top of a patch mesh

import os
import sys
import argparse
import shutil
from pathlib import Path

log_prefix = "[gen_dsl]"

def LOG(m):
    print("[gen_dsl]", m)

def ERROR(m):
    LOG(m)
    exit(1)


nadarepy_path = os.getenv('NADAREPY_PATH')
if nadarepy_path is None:
    LOG("Please set NADAREPY_PATH environment variable!"
        "(run 'source <build>/nadare_variables.sh')")
    exit(-1)
else:
    LOG("using nadarepy from " + nadarepy_path)

sys.path.insert(0, nadarepy_path)
import nadarepy

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=Path, help="input directory", required=True)
    parser.add_argument('-o', type=Path, help="output directory", required=True)
    parser.add_argument('-p', type=str, help="mesh patch name", required=True)
    parser.add_argument('-j', action="store_true", help="jitter points")
    parser.add_argument('--dt', type=float, help="time between samples in seconds")
    parser.add_argument('--duration', type=float, help="total duration in seconds")
    parser.add_argument('--cell-size', type=float, help="in meters")
    parser.add_argument('--velocity', type=float, help="in meters per second")
    # disk rotation parameters
    parser.add_argument('--rotating-disk', action="store_true", help="enable rotating disk")
    parser.add_argument('--disk-radius', type=float, help="in meters")
    parser.add_argument('--path-cx', type=float, help="in meters")
    parser.add_argument('--path-cy', type=float, help="in meters")
    parser.add_argument('--path-radius', type=float, help="in meters")
    # moving bar parameters
    parser.add_argument('--moving-bar', action="store_true", help="enable moving bar")
    parser.add_argument('--bar-cx', type=float, help="in meters")
    parser.add_argument('--bar-cy', type=float, help="in meters")
    parser.add_argument('--bar-width', type=float, help="in meters")
    parser.add_argument('--bar-depth', type=float, help="in meters")
    # stead front parameters
    parser.add_argument('--steady-front', action="store_true", help="enable steady front")
    parser.add_argument('--front-cx', type=float, help="in meters")
    parser.add_argument('--front-radius', type=float, help="in meters")
    args = parser.parse_args()

    if args.o.exists() and os.listdir(args.o):
        ERROR("[output folder][not empty] " + str(args.o))
    elif not args.o.exists():
        os.makedirs(args.o)

    LOG('-i: ' + str(args.i))
    LOG('-o: ' + str(args.o))
    LOG('-p: ' + str(args.p))                 
    LOG('-j:' + str(args.j) )
    LOG('--dt:' + str(args.dt))
    LOG('--duration: ' + str(args.duration))

    # simulation data
    of_sim = nadarepy.OFSim()
    of_sim.setSimPath(str(args.i))
    # setup procedural
    patch = of_sim.proceduralDSL(args.p)
    patch.setDuration(args.duration)
    patch.setWriteInterval(args.dt)
    # set rotating disk
    if args.rotating_disk:
        patch.addRotatingDisk(args.path_radius, args.disk_radius, args.velocity, args.path_cx, args.path_cy)
    elif args.moving_bar:
        patch.addMovingBar(args.bar_cx, args.bar_cy, args.bar_depth, args.bar_width, args.velocity)
    elif args.steady_front:
        patch.addSteadyFront(args.front_cx, args.front_radius, args.velocity)
    # output
    patch.write(str(args.o))

        
