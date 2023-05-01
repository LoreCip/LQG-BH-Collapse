import os
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



def GenerateVideo(ppath):
    X = np.loadtxt(ppath + '/outputs/xs.dat')
    B = np.loadtxt(ppath + '/outputs/B.dat')
    eps = np.loadtxt(ppath + '/outputs/E.dat')
    rho = np.loadtxt(ppath + '/outputs/rho.dat')
    times = np.loadtxt(ppath + '/outputs/times.dat')

    Path("./.frames/").mkdir(parents=True, exist_ok=True)

    def plotting(inputs):

        i, X, B, eps, rho, t = inputs 

        fig = plt.figure(figsize=(13, 10))
        fig.suptitle(f"Time: {np.round(t, 4)}", fontsize=14)
        ax1 = fig.add_subplot(3,1,1, adjustable='box')
        ax2 = fig.add_subplot(3,1,2)
        ax3 = fig.add_subplot(3,1,3)
        ax1.plot(X, rho)
        ax2.plot(X, B)
        ax3.plot(X, eps)

        ax1.set(xlim=[0,50], ylim=[0, 1.3*np.max(rho)])
        ax2.set(xlim=[0,50], ylim=[1.3*np.min(B), 0])
        ax3.set(xlim=[0,50], ylim=[1.3*np.min(eps), 1.3*np.max(eps)])

        ax2.set_ylabel("B")
        ax3.set_ylabel(r"$\epsilon^b$")
        ax1.set_ylabel(r"$\rho$")

        ax3.set_xlabel("x")

        plt.tight_layout()

        plt.savefig(f"./.frames/frame{i:04d}.png")
        plt.close()

    print('Creating frames.')
    inputs = [(i, X, B[i,:], eps[i,:], rho[i,:], times[i, 0]) for i in range(len(B))]
    p_map(plotting, inputs)


    print('Creating out_video.mp4')
    os.system("ffmpeg -framerate 30 -loglevel quiet -pattern_type glob -i '.frames/*.png' -c:v libx264 -pix_fmt yuv420p out_video.mp4")
    os.system("rm -rf ./.frames")


######################################################
###  FROM COMMAND LINE
######################################################

if __name__ == "__main__":

    usage = "usage: "+ sys.argv[0]+" -P [Path of simulation folder]"
    #################################
    # PARSE command line arguments
    #################################
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hP:")
    except :
        print(usage)
        sys.exit(2)

    for opt, arg in opts:
        if ( opt == '-h') :
            print(usage)
            sys.exit(2)
        elif (opt == '-P') :
            PATH = str(arg)

    ####################################
    # END PARSE command line arguments
    ####################################

    GenerateModel(PATH) 