"""This script is a general python script to plot results of 1D cases that are performed in 2D using D-Hydro FM"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import os
import configparser
from matplotlib.ticker import MaxNLocator, MultipleLocator
from matplotlib.widgets import Slider
import matplotlib as mpl
import math

# Set global font to Times New Roman, size 11, and mathtext italic
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.size'] = 11
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Times New Roman'
mpl.rcParams['mathtext.it'] = 'Times New Roman:italic'
mpl.rcParams['axes.titlesize'] = 11
mpl.rcParams['axes.labelsize'] = 11
mpl.rcParams['legend.fontsize'] = 11
mpl.rcParams['xtick.labelsize'] = 11
mpl.rcParams['ytick.labelsize'] = 11

# Set working directory

PATH = os.getcwd() # Ensure that this file is present in the test folder that contains ...\dflowfmoutput
OUTPUT_FOLDER = "dflowfmoutput" # The folder where the output files are stored
MDU_FILE_NAME = "kelvinfricno"
MAP_FILE_NAME = "kelvinfricno_map.nc"
FIGURES = PATH + '/Figures'
CASE_DIMENSION = 1 # Sets the case dimension. 1 = 1D -> All plots are f(x), 
                                            # 2 = 2D -> All plots are f(x,y) #Not implemented

# Figure settings
ms1 = 15    # markersize1
ms2 = 8     # markersize2
ms3 = 12    # markersize3
fs = 16     # fontsize
cl = 0.5    # color
lw = 1      # linewidth

# Load input parameters from mdu file
# Utility: Clean lines for configparser (strip comments, spaces)
def clean_mdu_lines(lines):
    cleaned = []
    for line in lines:
        # Remove comments after values
        if '#' in line:
            line = line.split('#', 1)[0]
        if '=' in line:
            key, value = line.split('=', 1)
            key = key.strip()
            value = value.strip()
            cleaned.append(f"{key} = {value}\n")
        else:
            cleaned.append(line)
    return cleaned

# Load and parse the mdu file
mdu_path = os.path.join(PATH, MDU_FILE_NAME + '.mdu')
config = configparser.ConfigParser()
config.optionxform = str  # preserve case

with open(mdu_path) as f:
    lines = f.readlines()
    lines = clean_mdu_lines(lines)
    config.read_string(''.join(lines))

# Convert configparser to nested dict (struct-like)
def config_to_struct(config):
    struct = {}
    for section in config.sections():
        struct[section] = {}
        for key, value in config.items(section):
            struct[section][key] = value
    return struct

mdu = config_to_struct(config)

# Example usage:
# mdu['geometry']['Uniformwidth1D']
# mdu['physics']['Ag']

# Optionally, convert values to float/int where possible
def try_cast(value):
    try:
        if '.' in value or 'e' in value.lower():
            return float(value)
        return int(value)
    except Exception:
        return value  # keep as string if not a number

for section in mdu:
    for key in mdu[section]:
        mdu[section][key] = try_cast(mdu[section][key])

class MDUSection:
    def __init__(self, entries):
        for key, value in entries.items():
            setattr(self, key, value)

class MDU:
    def __init__(self, mdu_dict):
        for section, entries in mdu_dict.items():
            setattr(self, section, MDUSection(entries))

# After you have built the mdu dictionary as before:
mdu_struct = MDU(mdu)

# Now you can access all entries as attributes, e.g.:
# mdu_struct.geometry.Uniformwidth1D
# mdu_struct.geometry.NetFile
# mdu_struct.model.Program

# uncomment below to print mdu entries for plotting implementation purposes
#for section in mdu:
#    print(f"[{section}]")
#    for key in mdu[section]:
#        print(f"  {key} = {getattr(getattr(mdu_struct, section), key)}")

dataset = nc.Dataset(PATH + '/' + OUTPUT_FOLDER + '/'+ MAP_FILE_NAME)

class NetCDFStruct:
    def __init__(self, dataset):
        # Load global attributes
        for attr in dataset.ncattrs():
            setattr(self, attr, getattr(dataset, attr))
        # Load variables (as netCDF4.Variable objects)
        for var in dataset.variables:
            setattr(self, var, dataset.variables[var])
        # Optionally, load dimensions
        self.dimensions = {dim: dataset.dimensions[dim].size for dim in dataset.dimensions}

# Usage:
nc_struct = NetCDFStruct(dataset)

# Print all variable names and shapes
print("Variables:")
for var in dataset.variables:
    print(f"  {var}: shape {dataset.variables[var].shape}")

# Print all global attributes
print("Global attributes:")
for attr in dataset.ncattrs():
    print(f"  {attr}: {getattr(dataset, attr)}")

# Load all data into structured arrays
# Automatically load all variables into a dictionary for easy access
data_vars = {}
for var in dataset.variables:
    data_vars[var] = dataset.variables[var][:]

# Plots

x = data_vars['mesh2d_face_x'][:]  # x-coordinates of faces
y = data_vars['mesh2d_face_y'][:]  # y-coordinates of faces
wd = data_vars['mesh2d_waterdepth'][:]  # water depth (time, face)
wl = data_vars['mesh2d_s1'][:] if 'mesh2d_s1' in data_vars else None  # water level (time, face)
vel = data_vars['mesh2d_ucmag'][:] if 'mesh2d_ucmag' in data_vars else None  # velocity magnitude (time, face)
times = data_vars['time'][:] if 'time' in data_vars else np.arange(wd.shape[0])

# --- Domain splitting logic based on y_max ---
y_max = y.max()
domain_bounds = [0, y_max/3, 2*y_max/3, y_max + 1e-6]  # Add small epsilon to include y_max itself
domain_indices = [
    np.where((y >= domain_bounds[i]) & (y < domain_bounds[i+1]))[0]
    for i in range(3)
]

# --- Analytical solution for water level ζ(y, t) ---
# Model parameters from the screenshot
B = 550e3  # width in meters
H = 80     # water depth in meters
g = 9.81   # gravitational acceleration in m/s^2
phi = 52 * np.pi / 180  # latitude in radians
Omega = 2 * np.pi / (24*3600)  # Earth's rotation rate in rad/s
f = 2 * Omega * np.sin(phi)    # Coriolis parameter
c = np.sqrt(g * H)             # phase speed
U0 = c / g                     # forcing amplitude
T = 745 * 60                   # wave period in seconds
omega = 2 * np.pi / T          # wave frequency
k = omega / c                  # wavenumber

def analytical_zeta(x, y, t):
    """
    Analytical solution for water level ζ(x, y, t) as a superposition of two Kelvin waves,
    following Roos et al. (2011), eq. (A.7), with boundary forcing at x = L.
    x, y in meters, t in seconds.
    """
    L = x.max()  # or set L to your known domain length
    A = (2 * c / g) * U0 * np.exp(-f * B / (2 * c))
    exp_fac = np.exp(-f * y / c)
    zeta = A * exp_fac * (np.cos(k * (L - x) - omega * t) + np.cos(k * (L + x) - omega * t))
    return zeta

# --- Plotting as contour plots ---
fig, axes = plt.subplots(4, 2, figsize=(11, 8))
plt.subplots_adjust(bottom=0.15, hspace=0.35)

contour_sets = []
colorbars = []
time_idx = 0  # Initial time index

# Domain names in requested order: domain 3, 2, 1
domain_names = ['triangle', 'squares_fine', 'squares']

for plot_idx, domain_idx in enumerate([2, 1, 0]):
    inds = domain_indices[domain_idx]
    # Water level contour
    ax_wl = axes[plot_idx, 0]
    if wl is not None and len(inds) > 0:
        c_wl = ax_wl.tricontourf(
            x[inds], y[inds], wl[time_idx, inds], levels=20, cmap='Blues'
        )
        cb_wl = fig.colorbar(c_wl, ax=ax_wl, orientation='vertical', pad=0.01)
        cb_wl.set_label('Water Level')
        colorbars.append(cb_wl)
        ax_wl.set_ylabel(f'{domain_names[plot_idx]}\nWater Level')
        contour_sets.append(c_wl)
    else:
        ax_wl.text(0.5, 0.5, 'No water level data', ha='center', va='center')
        contour_sets.append(None)
    ax_wl.set_aspect('equal')
    ax_wl.grid(True)

    # Velocity magnitude contour
    ax_vel = axes[plot_idx, 1]
    if vel is not None and len(inds) > 0:
        c_vel = ax_vel.tricontourf(
            x[inds], y[inds], vel[time_idx, inds], levels=20, cmap='Reds'
        )
        cb_vel = fig.colorbar(c_vel, ax=ax_vel, orientation='vertical', pad=0.01)
        cb_vel.set_label('Velocity')
        colorbars.append(cb_vel)
        ax_vel.set_ylabel(f'{domain_names[plot_idx]}\nVelocity')
        contour_sets.append(c_vel)
    else:
        ax_vel.text(0.5, 0.5, 'No velocity data', ha='center', va='center')
        contour_sets.append(None)
    ax_vel.set_aspect('equal')
    ax_vel.grid(True)

# --- Analytical solution row (4th row, left column only) ---
ax_analytical = axes[3, 0]
analytical_cbar = None  # Store the colorbar for the analytical plot
# Only plot for the fine rectangular domain (domain_idx = 1)
inds = domain_indices[1]
if len(inds) > 0:
    x_domain = x[inds]
    y_domain = y[inds]
    zeta_analytical = analytical_zeta(x_domain, y_domain, times[time_idx])
    c_analytical = ax_analytical.tricontourf(
        x_domain, y_domain, zeta_analytical, levels=20, cmap='Blues'
    )
    analytical_cbar = fig.colorbar(c_analytical, ax=ax_analytical, orientation='vertical', pad=0.01)
    analytical_cbar.set_label('Analytical Water Level')
    ax_analytical.set_ylabel('Analytical\nWater Level')
    ax_analytical.set_xlabel('x / y')
    ax_analytical.set_title('Analytical Water Level (Fine Rectangular Domain)')
    ax_analytical.set_aspect('equal')
    ax_analytical.grid(True)
else:
    ax_analytical.text(0.5, 0.5, 'No analytical data', ha='center', va='center')

# Leave the 4th row, 2nd column empty
axes[3, 1].axis('off')

axes[-2, 0].set_xlabel('x')
axes[-2, 1].set_xlabel('x')
axes[0, 0].set_title('Water Level')
axes[0, 1].set_title('Velocity Magnitude')
axes[3, 0].set_title('Analytical Water Level')

# --- Slider ---
ax_slider = plt.axes([0.2, 0.05, 0.6, 0.03])
slider = Slider(ax_slider, 'Time', 0, len(times)-1, valinit=0, valstep=1)

def update(val):
    global analytical_cbar
    idx = int(slider.val)
    # Remove all colorbars before redrawing
    for cb in colorbars:
        cb.remove()
    colorbars.clear()
    if analytical_cbar is not None:
        analytical_cbar.remove()
        analytical_cbar = None
    for plot_idx, domain_idx in enumerate([2, 1, 0]):
        inds = domain_indices[domain_idx]
        # Update water level contour
        ax_wl = axes[plot_idx, 0]
        for coll in ax_wl.collections:
            coll.remove()
        if wl is not None and len(inds) > 0:
            c_wl = ax_wl.tricontourf(
                x[inds], y[inds], wl[idx, inds], levels=20, cmap='Blues'
            )
            cb_wl = fig.colorbar(c_wl, ax=ax_wl, orientation='vertical', pad=0.01)
            cb_wl.set_label('Water Level')
            colorbars.append(cb_wl)
            contour_sets[2*plot_idx] = c_wl
        # Update velocity contour
        ax_vel = axes[plot_idx, 1]
        for coll in ax_vel.collections:
            coll.remove()
        if vel is not None and len(inds) > 0:
            c_vel = ax_vel.tricontourf(
                x[inds], y[inds], vel[idx, inds], levels=20, cmap='Reds'
            )
            cb_vel = fig.colorbar(c_vel, ax=ax_vel, orientation='vertical', pad=0.01)
            cb_vel.set_label('Velocity')
            colorbars.append(cb_vel)
            contour_sets[2*plot_idx+1] = c_vel
    # Update analytical plot
    ax_analytical = axes[3, 0]
    ax_analytical.cla()
    inds = domain_indices[1]
    if len(inds) > 0:
        x_domain = x[inds]
        y_domain = y[inds]
        zeta_analytical = analytical_zeta(x_domain, y_domain, times[idx])
        c_analytical = ax_analytical.tricontourf(
            x_domain, y_domain, zeta_analytical, levels=20, cmap='Blues'
        )
        analytical_cbar = fig.colorbar(c_analytical, ax=ax_analytical, orientation='vertical', pad=0.01)
        analytical_cbar.set_label('Analytical Water Level')
        ax_analytical.set_ylabel('Analytical\nWater Level')
        ax_analytical.set_xlabel('x / y')
        ax_analytical.set_title('Analytical Water Level (Fine Rectangular Domain)')
        ax_analytical.set_aspect('equal')
        ax_analytical.grid(True)
    else:
        ax_analytical.text(0.5, 0.5, 'No analytical data', ha='center', va='center')
    axes[3, 1].axis('off')
    fig.suptitle(f'Time step: {idx} (t={times[idx]:.2f})', fontsize=14)
    fig.canvas.draw_idle()

slider.on_changed(update)
fig.suptitle(f'Time step: {time_idx} (t={times[time_idx]:.2f})', fontsize=14)
plt.show()
