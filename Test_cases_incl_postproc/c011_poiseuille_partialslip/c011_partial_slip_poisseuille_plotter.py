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
MDU_FILE_NAME = "poiseuillepartialslip"
MAP_FILE_NAME = "poiseuillepartialslip_map.nc"
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
#wd = data_vars['mesh2d_waterdepth'][:]  # water depth (time, face)
wl = data_vars['mesh2d_s1'][:] if 'mesh2d_s1' in data_vars else None  # water level (time, face)
vel = data_vars['mesh2d_ucmag'][:] if 'mesh2d_ucmag' in data_vars else None  # velocity magnitude (time, face)
# times = data_vars['time'][:] if 'time' in data_vars else np.arange(wd.shape[0])

import numpy as np

# --- Parameters from Thacker test case ---
D0 = 10.0        # reference depth [m]
eta = 2.0        # max surface elevation [m]
L = 102000.0     # paraboloid radius [m]
g = 9.81         # gravity [m/s^2]

# Compute A and omega
A = ((D0 + eta)**2 - D0**2) / ((D0 + eta)**2 + D0**2)
omega = np.sqrt(g * D0) / L

# --- Analytical solution functions ---

def analytical_velocity(g, dH, Dy, B, L, alpha, nu_h):
    c0 = -g*dH/(2*L)
    c1 = nu_h/alpha * sqrt(B*c0)+ c0*((B*Dy)/2)**2
    return (c1-c0*y**2)/nu_h


# Add this setting at the top of your plotting section:
PLOT_COMBINED = True  # Set to False to plot individually, True for combined plot
PLOT_TEMPORAL_EVOLUTION = False  # Set to True to enable the temporal evolution plot

from matplotlib.widgets import Slider

    
# --- Parameters ---
nu_h = 0.1
DH = -1e-4
L = 10000.0   # [m]
B = 1000.0    # [m]
kappa = 0.41

# Get k_s from MDU file (assume it's in mdu_struct.physics.ks)
k_s = getattr(mdu_struct.physics, 'wall_ks', 0.1)
y_0 = k_s / 30

# --- Select cross-section at x = midpoint ---
x_section = np.mean(x)  # or set to desired x value
cross_idx = np.where(x == 9925.0)[0]
x_cross = x[cross_idx]
y_cross = y[cross_idx]
u_num = data_vars['mesh2d_ucx'][-1, cross_idx]  # t=0, adjust index as needed

# --- Analytical solution ---
Dy = np.mean(np.diff(np.sort(y_cross)))  # grid size along cross-section
alpha = kappa / np.log(1 + Dy / (2 * y_0))

def analytical_u_velocity(y, nu_h, alpha, DH, B, L, g=9.81):
    # y: position across channel (centered at 0)
    c0 = -g * DH / (2 * L)
    c1 = nu_h / alpha * np.sqrt(B * c0) + c0 * ((B - Dy) / 2) ** 2
    return (c1 - c0 * y ** 2) / nu_h

# For a cross-section, y goes from -B/2 to B/2
y_analytical = y_cross - np.mean(y_cross)  # center at 0
u_analytical = analytical_u_velocity(y_analytical, nu_h, alpha, DH, B, L)
a = 1 + 1

# --- Plot ---
import matplotlib.pyplot as plt

plt.figure(figsize=(5, 5))
plt.plot(u_num, y_cross-B/2, 'bo-', label='Numerical u-velocity')
plt.plot(u_analytical, y_cross-B/2, 'r--', label='Analytical u-velocity')
plt.xlabel('u-velocity [m/s]')
plt.ylabel('y [m] (cross-section)')
plt.title('Cross-section u-velocity comparison')
plt.legend()
plt.grid(True, linestyle=':')
plt.tight_layout()
plt.show()

relative_error = (np.max(u_analytical) - np.max(u_num)) / np.max(u_analytical)
print(f"Relative error: {relative_error:.2%}")
