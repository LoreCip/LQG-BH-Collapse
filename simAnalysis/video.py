import sys
import os
import warnings

import h5py
import readSim as rs

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

def determine_y_axis_range(data):
    min_value = np.min(data)
    max_value = np.max(data)

    if min_value > 0:
        y_min = 0
        y_max = max_value * 1.3  # Add some padding above the maximum value
    elif max_value < 0:
        y_min = min_value * 1.3  # Add some padding below the minimum value
        y_max = 0
    else:
        y_min = min_value * 1.3  # Add some padding below the minimum value
        y_max = max_value * 1.3  # Add some padding above the maximum value

    return [y_min, y_max]

def plotting(inputs):

        i, X, B, eps, rho, t, m = inputs 

        fig = plt.figure(figsize=(13, 10))
        fig.suptitle(f'Time: {np.round(t, 4)}', fontsize=14)
        ax1 = fig.add_subplot(3,1,1, adjustable='box')
        ax2 = fig.add_subplot(3,1,2)
        ax3 = fig.add_subplot(3,1,3)
        ax1.plot(X, rho)
        ax2.plot(X, B)
        ax3.plot(X, eps)

        xlim=[0, 2*m*2.5]

        ax1.set(xlim=xlim, ylim=determine_y_axis_range(rho[100:-100]))
        ax2.set(xlim=xlim, ylim=determine_y_axis_range(B[100:-100]))
        ax3.set(xlim=xlim, ylim=determine_y_axis_range(eps[100:-100]))

        ax2.set_ylabel('B')
        ax3.set_ylabel(r'$\epsilon^b$')
        ax1.set_ylabel(r'$\rho$')

        ax3.set_xlabel('x')

        plt.tight_layout()

        plt.savefig(f'./.frames/frame{i:07d}.png')
        plt.close()

def ProduceInputs(sim, last_line, pline):

    it = np.array(range(last_line, last_line + pline))
    it = it[it <= sim.niter]

    X = sim.xgrid
    Bs = sim.get(it, 'B')
    ebs = sim.get(it, 'e^b')
    rhos = sim.get(it, 'rho')
    ts = sim.get(it, 't')
    m = sim.mass

    return [(last_line + i, X, Bs[i], ebs[i], rhos[i], ts[i], m) for i in range(last_line, last_line + pline)]

def GenerateVideo(ppath, n_partitions):

    if n_partitions <= 0:
        n_partitions = 1

    Path("./.frames/").mkdir(parents=True, exist_ok=True)

    sim = rs.Sim(ppath)

    print('Creating frames.')
    
    last_line = 0
    p_line = int( sim.niter / n_partitions)
    
    for j in range(n_partitions):
        inputs = ProduceInputs(sim, last_line, p_line)
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