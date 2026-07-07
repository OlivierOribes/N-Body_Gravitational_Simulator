"""
@file animation_2D_innermost_planet_leapfrog.py
@author Olivier Oribes
@brief Generate a 2D animation comparing simulated and NASA Horizons planetary trajectories.

This script reads Cartesian position files produced by the N-body Solar System
simulation and generates a 2D MP4 animation of the selected planets.

Two trajectory datasets are displayed simultaneously:
  - Numerical trajectories computed by the simulation
  - Reference trajectories extracted from NASA JPL Horizons

The calculated trajectories are drawn with solid lines, while NASA reference
trajectories are drawn with dashed lines. Planet radii are loaded from
settingdata.txt to scale the marker sizes approximately for visualization.

Input files:
  - settingdata.txt         : body names and radii
  - positions_leapfrog.txt  : simulated Cartesian positions
  - positionsnasa.txt       : NASA Horizons reference positions

Output file:
  - animation_2D_innermost_planet_leapfrog.mp4

Main features:
  - 2D orbital animation in the XY plane
  - Comparison between simulated and NASA reference trajectories
  - Separate color maps for numerical and NASA data
  - Planet filtering through allowed_planets
  - MP4 export using Matplotlib and FFmpeg

Target system:
  - Python 3
  - Matplotlib
  - FFmpeg installed and available in PATH

@date 2026-07-06
@copyright Copyright (c) 2026 Olivier Oribes
@license BSD-3-Clause
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from collections import defaultdict

# ==============================================================
#  Define Color of planets (with a class dictionary)
# ==============================================================
# Colors for bodies from calculated positions (with solid lines)

calc_color_map = {
    "Sun":      "yellow",
    "Mercury":  "mediumvioletred",
    "Venus":    "gray",
    "Earth":    "blue",
    "Mars":     "red",
    "Jupiter":  "brown",
    "Saturn":   "gold",
    "Uranus":   "cyan",
    "Neptune":  "darkblue"
}

# Colors for bodies from NASA positions (with dashed lines)
nasa_color_map = {
    "Sun":      "goldenrod",
    "Mercury":  "purple",
    "Venus":    "black",
    "Earth":    "green",
    "Mars":     "orange",
    "Jupiter":  "maroon",
    "Saturn":   "navy",
    "Uranus":   "teal",
    "Neptune":  "indigo"
}

# ==============================================================
#  Define Allowed Planets to Plot
# ==============================================================
# Only these planets will be plotted (use to reduce time simulation
# or to focus on a particular orbital)
allowed_planets = {"Sun", "Mercury", "Venus", "Earth", "Mars"} 

# ==============================================================
# Load Planet Sizes from settingdata.txt
# ==============================================================

# Could be disabled but it shows planets' scale difference
planet_sizes = {}

with open('settingdata.txt', 'r') as f:

    for line in f:
        line = line.strip()

        if line.startswith('#') or not line:
            continue

        parts = line.split()

        if len(parts) == 2:
            name, size = parts
            planet_sizes[name] = float(size)

# ==============================================================
# Function to Load Positions from a File
# ==============================================================
def load_positions(filename):

    # Reads a positions file and returns two objects:
    # - data: a dictionary mapping each simulation step (int) to a list of tuples (name, x, y, z)
    # - steps: a sorted list of simulation step numbers (that represent days)
    #
    # Needs the file to start with "# Step <number>" followed by one line per body
    # in the format:
    #     << Name   x   y   z >>

    data = {}
    steps = []
    with open(filename, 'r') as f:
        current_step = None
        for line in f:
            line = line.strip()

            if line.startswith("# Step"):
                parts = line.split()
                current_step = int(parts[2])
                steps.append(current_step)
                data[current_step] = []
            elif line:
                parts = line.split()
                if len(parts) >= 4:
                    name = parts[0]
                    x, y, z = map(float, parts[1:4])
                    data[current_step].append((name, x, y, z))
    return data, sorted(steps)

# ==============================================================
# Load Position Data from Files
# ==============================================================
# Load calculated positions (often something like "positions.txt")
calc_data, calc_steps = load_positions('positions_leapfrog.txt')
# Load official NASA positions
nasa_data, nasa_steps = load_positions('positionsnasa.txt')
total_steps = len(calc_steps)
steps = calc_steps

# ==============================================================
# Build Trajectories for Each Set
# ==============================================================
# For calculated positions: build a dictionary mapping planet :
#     name -> (list of xs, list of ys, list of zs)
calc_trajectories = defaultdict(lambda: ([], [], []))
for step in steps:
    for name, x, y, z in calc_data[step]:
        if name not in allowed_planets: # reduce simulation time by restricting to allowed planets
            continue
        calc_trajectories[name][0].append(x)
        calc_trajectories[name][1].append(y)
        calc_trajectories[name][2].append(z)

# For NASA positions: same as above
nasa_trajectories = defaultdict(lambda: ([], [], []))
for step in nasa_steps:
    for name, x, y, z in nasa_data[step]:
        if name not in allowed_planets:
            continue
        nasa_trajectories[name][0].append(x)
        nasa_trajectories[name][1].append(y)
        nasa_trajectories[name][2].append(z)

# ==============================================================
# Set Up the Plot
# ==============================================================
fig, ax = plt.subplots()
max_range = 1.8
ax.set_xlim([-max_range, max_range])
ax.set_ylim([-max_range, max_range])
ax.set_aspect('equal', adjustable='box')

# ==============================================================
# Update Function for the Animation
# ==============================================================

def update(frame):
    ax.clear()
    ax.set_xlim([-max_range, max_range])
    ax.set_ylim([-max_range, max_range])
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(f"Step {frame}")

    for planet, (xs, ys, _) in calc_trajectories.items():

        # Sun is sized separately: its true scale would dwarf the other planets
        if planet == "Sun":
            size = max(planet_sizes.get(planet, 1) * 0.0001, 1)
        else:
            size = max(planet_sizes.get(planet, 1) * 0.02, 2)

        if planet not in allowed_planets:
            continue
        color_calc = calc_color_map.get(planet, "white")
        if frame < len(xs):
            label = f"{planet} calc"
            ax.plot(xs[:frame+1], ys[:frame+1], color=color_calc, linewidth=0.5, label=label)
            # Redraws the whole trajectory every frame, not just incrementally - not efficient.
            ax.scatter(xs[frame], ys[frame], s=size, color=color_calc)

    for planet, (xs, ys, _) in nasa_trajectories.items():

        if planet == "Sun":
            size = max(planet_sizes.get(planet, 1) * 0.0001, 1)
        else:
            size = max(planet_sizes.get(planet, 1) * 0.02, 2)

        if planet not in allowed_planets:
            continue
        color_nasa = nasa_color_map.get(planet, "white")
        if frame < len(xs):
            label = f"{planet} NASA"
            ax.plot(xs[:frame+1], ys[:frame+1], color=color_nasa, linestyle="--", linewidth=0.5, label=label)
            ax.scatter(xs[frame], ys[frame], s=size, color=color_nasa)

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        unique = dict(zip(labels, handles)) # dict dedupes labels while pairing them with their handle
        sorted_labels = sorted(unique.keys())
        sorted_handles = [unique[label] for label in sorted_labels]

        ax.legend(sorted_handles, sorted_labels, fontsize=8, loc='upper left',
                  bbox_to_anchor=(1.05, 1), borderaxespad=0.)


# ==============================================================
# Create the Animation
# ==============================================================
frames_to_plot = range(0, 63100, 5)
ani = FuncAnimation(fig, update, frames=frames_to_plot, interval=50, blit=False)
# ==============================================================
# Save the Animation as MP4
# ==============================================================
writer = FFMpegWriter(fps=10, bitrate=1800)
ani.save("animation_2D_innermost_planet_leapfrog.mp4", writer=writer)

print("Animation saved as 'animation_2D_innermost_planet_leapfrog.mp4'")
