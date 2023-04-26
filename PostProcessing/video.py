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

X = np.loadtxt('/home/lorenzo/phd/LQG/negative_build/outputs/xs.dat')
B = np.loadtxt('/home/lorenzo/phd/LQG/negative_build/outputs/B.dat')
eps = np.loadtxt('/home/lorenzo/phd/LQG/negative_build/outputs/E.dat')
rho = np.loadtxt('/home/lorenzo/phd/LQG/negative_build/outputs/rho.dat')


Path("./.frames/").mkdir(parents=True, exist_ok=True)

for i in range(len(B)):
    fig = plt.figure()
    plt.plot(X[0,:], B[i, :])
    plt.savefig(f"./.frames/frame{i:04d}.png")
    plt.close()

print('Creating out_video.mp4')
os.system("ffmpeg -framerate 30 -loglevel quiet -pattern_type glob -i '.frames/*.png' -c:v libx264 -pix_fmt yuv420p out_video.mp4")
os.system("rm -rf ./.frames")