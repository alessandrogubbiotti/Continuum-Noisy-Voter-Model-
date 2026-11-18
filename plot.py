"""plot.py -- draws the colours written by nvm_sim.c.

Reads output.txt (a first line with the parameters, then one line for each
saved time with the colours of the sites) and saves output.png: space on the
horizontal axis, time on the vertical one, the two colours of the sites.

    python3 plot.py
"""
import matplotlib
matplotlib.use("Agg")          # no window is needed: the picture is saved to a file
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

# The first line of the file carries the parameters of the run.
with open("output.txt") as f:
    header = f.readline().strip().lstrip("# ")
parameters = dict(item.split("=") for item in header.split())
colours = np.loadtxt("output.txt", ndmin=2)

# The row j is the time j/resolution: it is drawn as a band of height
# 1/resolution, centred at its time.
length = float(parameters["L"])
step = 1.0 / float(parameters["resolution"])
last_time = (colours.shape[0] - 1) * step

plt.figure(figsize=(10, 4))
plt.imshow(colours, cmap=ListedColormap(["#2166ac", "#d73027"]), vmin=0, vmax=1,
           aspect="auto", interpolation="nearest", origin="lower",
           extent=(0, length, -step / 2, last_time + step / 2))
plt.ylim(0, max(last_time, step / 2))
plt.xlabel("Space")
plt.ylabel("Time")
plt.title("Noisy Voter Model  (" + header + ")")

plt.tight_layout()
plt.savefig("output.png", dpi=150)
