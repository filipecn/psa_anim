# this script processes all statistics data from simulations in a given directory
# and generates a report
# the stats file is composed of a time series containing data
# each line of a stats file contains:
# 0 time 
# 1 mass
# 2 entrainment
# 3 min u
# 4 max u
# 5 kinetic energy 
# 6 energy entrainment

import os
import argparse
import shutil
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as mp

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=Path, help="sim directory")
    parser.add_argument("-o", type=Path, help="report output directory")
    parser.add_argument("-n", type=str, help="stats file base name")
    parser.add_argument('-s', action='store_true', help="display figures")
    parser.add_argument('-w', type=int, help='figure width in pixels', default=1024)
    parser.add_argument('--height', type=int, help='figure height in pixels', default=1024)
    parser.add_argument('--dpi', type=int, help='monitor dpi', default=96)
    args = parser.parse_args()
    
    # group experiments by type (lock_exchange, entrainment, full)
    groups = {}

    # go over all directories and find stats files
    for dirpath, dirs, files in os.walk(args.i):
        sim_type = os.path.basename(dirpath)
        # check for stats file
        if args.n + "0" in files:
            if sim_type not in groups:
                groups[sim_type] = [dirpath]
            else:
                groups[sim_type].append(dirpath)

    # for each group, generate graphs containing curves from all simulations in the group
    for group_name in groups:
        # for each simulation path retrieve the data from the stats files
        group_data = []
        group_sim_names = []
        for sim_path in groups[group_name]:
            print(sim_path)
            # find out how many stat files there are (parallel runs produce multiple files)
            stats_files = [f for f in os.listdir(sim_path) if args.n in f and 'labels' not in f]
            # get data name
            labels = []

            print('/'.join([sim_path, args.n + '_labels']))
            with open('/'.join([sim_path, args.n + '_labels']), 'r') as f:
                labels = f.readline().split()
            # read data 
            df = pd.DataFrame()
            for stats_file in stats_files:
                print('/'.join([sim_path, stats_file])) 
                data = pd.read_csv('/'.join([sim_path, stats_file]), sep=' ', 
                        names=labels, header=None, index_col=False)
                print(data)
                if len(df) == 0:
                    df = data
                else:
                    for c in df.columns[1:len(df.columns)+1]:
                        df.loc[:,c] = [a+b for a, b in zip(df.loc[:,c],data.loc[:,c])]
            # register group
            group_data.append(df)
            group_sim_names.append(sim_path.split('/')[-4])
        
        # gen plots
        #
        #   mass        mass entr
        #
        #   energy      energy entr
        #
        fig, axes = mp.subplots(ncols=2, nrows=2, figsize=(args.w / args.dpi, args.height / args.dpi))
        ###########################################################################################
        # TOTAL MASS                                                                              #
        ###########################################################################################
        axx = axes[0,0]
        for i in range(len(group_data)):
            group_data[i].plot(ax=axx, x=0, y=['total_mass'], kind='line', label=[group_sim_names[i]])
        axx.set_xlabel('time (s)')
        axx.set_ylabel('total mass (kg)')
        ###########################################################################################
        # MASS ENTRAINMENT                                                                        #
        ###########################################################################################
        axx = axes[0,1]
        for i in range(len(group_data)):
            group_data[i].plot(ax=axx, x=0, y=['total_entrained_mass'], kind='line', 
                    label=[group_sim_names[i]])
        for i in range(len(group_data)):
            # here we plot the diff between subsequent time steps
            ddf = group_data[i][['total_mass']].diff()
            ddf['time'] = group_data[i]['time']
            ddf.plot(ax=axx, x='time', y=['total_mass'], style=['--'], kind='line', 
                    label=['D ' + group_sim_names[i]])
        axx.set_xlabel('time (s)')
        axx.set_ylabel('entrained mass (kg)')
        ###########################################################################################
        # TOTAL ENERGY                                                                            #
        ###########################################################################################
        axx = axes[1,0]
        for i in range(len(group_data)):
            group_data[i].plot(ax=axx, x=0, y=['total_kinetic_energy'], kind='line', 
                    label=['K ' + group_sim_names[i]])
        for i in range(len(group_data)):
            group_data[i].plot(ax=axx, x=0, y=['total_potential_energy'], style=['--'], kind='line', 
                    label=['P ' + group_sim_names[i]])
        axx.set_xlabel('time (s)')
        axx.set_ylabel('energy (kg m2 /s2)')
        ###########################################################################################
        # ENERGY ENTRAINMENT                                                                      #
        ###########################################################################################
        axx = axes[1,1]
        for i in range(len(group_data)):
            group_data[i].plot(ax=axx, x=0, y=['total_entrained_kinetic_energy'], kind='line', 
                    label=[group_sim_names[i]])
        for i in range(len(group_data)):
            # here we plot the diff between subsequent time steps
            ddf = group_data[i][['total_kinetic_energy']].diff()
            ddf['time'] = group_data[i]['time']
            ddf.plot(ax=axx, x='time', y=['total_kinetic_energy'], style=['--'], kind='line', 
                    label=['D ' + group_sim_names[i]])
        axx.set_xlabel('time (s)')
        axx.set_ylabel('energy (kg m2 /s2)')
        fig.suptitle(group_name)
        
        # store graph
        print(args.o / (group_name + '.png'))
        mp.savefig(args.o / (group_name + '.png'))
        if args.s:
            mp.show()
                    
