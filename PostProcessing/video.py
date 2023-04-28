import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

########
## TO DO
########
# 1) check for existance of frame folder and change name
# 2) check for existance of out_video and change name
# 3) input name of out_video
# 4) option to keep frames

X = np.loadtxt('/home/lorenzo/phd/LQG/nbuild/outputs/xs.dat')
B = np.loadtxt('/home/lorenzo/phd/LQG/nbuild/outputs/B.dat')
eps = np.loadtxt('/home/lorenzo/phd/LQG/nbuild/outputs/E.dat')
rho = np.loadtxt('/home/lorenzo/phd/LQG/nbuild/outputs/rho.dat')


Path("./.frames/").mkdir(parents=True, exist_ok=True)

for i in range(len(B)):
    fig = plt.figure(figsize=(13, 10))
    ax1 = fig.add_subplot(3,1,1, adjustable='box')
    ax2 = fig.add_subplot(3,1,2)
    ax3 = fig.add_subplot(3,1,3)
    ax1.plot(X[:], rho[i, :])
    ax2.plot(X[:], B[i, :])
    ax3.plot(X[:], eps[i, :])

    ax1.set(xlim=[0,50], ylim=[-1e-4, 1.3*np.max(rho[i,:])])
    ax2.set(xlim=[0,50], ylim=[1.3*np.min(B[i, :]), 0])
    ax3.set(xlim=[0,50], ylim=[0, 1.3*np.max(eps[i,:])])

    ax2.set_ylabel("B")
    ax3.set_ylabel(r"$\epsilon^b$")
    ax1.set_ylabel(r"$\rho$")

    ax3.set_xlabel("x")

    plt.tight_layout()

    plt.savefig(f"./.frames/frame{i:04d}.png")
    plt.close()

print('Creating out_video.mp4')
os.system("ffmpeg -framerate 15 -loglevel quiet -pattern_type glob -i '.frames/*.png' -c:v libx264 -pix_fmt yuv420p out_video.mp4")
os.system("rm -rf ./.frames")