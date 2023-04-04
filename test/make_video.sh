#!/bin/bash

ffmpeg -framerate 30 -loglevel quiet -pattern_type glob -i 'frames/*.png' -c:v libx264 -pix_fmt yuv420p out_video.mp4