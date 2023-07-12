import sys
import os
import warnings

import h5py

import getopt
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from p_tqdm import p_map

########
## TO DO
########
# 1) check for existance of frame folder and change name
# 2) check for existance of out_video and change name
# 3) input name of out_video
# 4) option to keep frames

def plotting(inputs):

        i, X, B, eps, rho, t = inputs 

        fig = plt.figure(figsize=(13, 10))
        fig.suptitle(f'Time: {np.round(t, 4)}', fontsize=14)
        ax1 = fig.add_subplot(3,1,1, adjustable='box')
        ax2 = fig.add_subplot(3,1,2)
        ax3 = fig.add_subplot(3,1,3)
        ax1.plot(X, rho)
        ax2.plot(X, B)
        ax3.plot(X, eps)

        ax1.set(xlim=[0,50], ylim=[0, min(3, 1.3*np.max(rho))])
        ax2.set(xlim=[0,50], ylim=[1.3*np.min(B), 0])
        ax3.set(xlim=[0,50], ylim=[1e-2+1.3*np.min(eps), 1.3*np.max(eps)])

        ax2.set_ylabel('B')
        ax3.set_ylabel(r'$\epsilon^b$')
        ax1.set_ylabel(r'$\rho$')

        ax3.set_xlabel('x')

        plt.tight_layout()

        plt.savefig(f'./.frames/frame{i:07d}.png')
        plt.close()

def count_h5_groups(file_path):
    count = 0

    with h5py.File(file_path, 'r') as f:
        # Recursive function to traverse the HDF5 file structure
        def traverse_groups(group):
            nonlocal count
            count += 1  # Increment count for each group encountered

            for key in group.keys():
                if isinstance(group[key], h5py.Group):
                    traverse_groups(group[key])  # Recursively visit subgroups

        # Start traversing from the root group
        traverse_groups(f)

    return count

def ProduceInputs(ppath, names, last_line, pline):

    with h5py.File(ppath+'/outputs/output.h5', 'r') as f:
        X = f['Xgrid'][()]

        if len(names) == 1:
            out = [(last_line + 1, X, f[str(names[0])]['B'][()], f[str(names[0])]['e^b'][()], f[str(names[0])]['rho'][()], f[str(names[0])]['t'][()])]
        else:
            out = [(last_line + i, X, f[str(names[i])]['B'][()], f[str(names[i])]['e^b'][()], f[str(names[i])]['rho'][()], f[str(names[i])]['t'][()]) for i in range(len(names))]
        
    return out

def GenerateVideo(ppath, n_partitions):

    if n_partitions <= 0:
        n_partitions = 1

    Path("./.frames/").mkdir(parents=True, exist_ok=True)

    print('Creating frames.')
    
    last_line = 0
    p_line = int( count_h5_groups(ppath+'/outputs/output.h5') / n_partitions)

    with h5py.File(ppath+'/outputs/output.h5', 'r') as f:
        group_names = list(f.keys())
        rl = []
        for i in range(len(group_names)):
            try:
                int(group_names[i])  # Try converting the element to an integer
            except ValueError:
                rl.append(i)  # Remove the element if it cannot be converted
        for r in rl[::-1]:
            group_names.pop(r)

    group_names = sorted(group_names, key=int)

    for j in range(n_partitions):
        inputs = ProduceInputs(ppath, group_names[last_line:last_line+p_line], last_line, p_line)
        p_map(plotting, inputs)
        last_line += p_line

    print('Creating out_video.mp4')
    os.system("ffmpeg -framerate 30 -loglevel quiet -pattern_type glob -i '.frames/*.png' -c:v libx264 -pix_fmt yuv420p out_video.mp4")
    os.system("rm -rf ./.frames")


######################################################
###  FROM COMMAND LINE
######################################################

if __name__ == "__main__":

    batches = 1

    usage = "usage: "+ sys.argv[0]+" -P [Path of simulation folder] -b [Number of batches]"
    #################################
    # PARSE command line arguments
    #################################
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hP:b:")
    except :
        print(usage)
        sys.exit(2)

    for opt, arg in opts:
        if ( opt == '-h') :
            print(usage)
            sys.exit(2)
        elif (opt == '-P') :
            PATH = str(arg)
        elif (opt == '-b') :
            batches = int(arg)

    ####################################
    # END PARSE command line arguments
    ####################################

    GenerateVideo(PATH, batches) 