"""This script is a general python script to plot a cross-section of water depth at the middle of the domain along the width of the channel and compare it with the analytical Coriolis solution."""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import os

# Set global font to Times New Roman, size 11, and mathtext italic
import matplotlib as mpl
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
PATH = os.getcwd()
OUTPUT_FOLDER = "dflowfmoutput"
MAP_FILE_NAME = "flatbed_quadrilateral_map.nc"

# Load NetCDF data
dataset = nc.Dataset(os.path.join(PATH, OUTPUT_FOLDER, MAP_FILE_NAME))
data_vars = {var: dataset.variables[var][:] for var in dataset.variables}

x = data_vars['mesh2d_face_x'][:]
y = data_vars['mesh2d_face_y'][:]
wd = data_vars['mesh2d_waterdepth'][:]  # shape: (time, face)
u = data_vars['mesh2d_ucx'][:] if 'mesh2d_ucx' in data_vars else None
times = data_vars['time'][:] if 'time' in data_vars else np.arange(wd.shape[0])

# Time setup
t_index = 5
t = times[t_index]

# --- Cross-section plot at middle of domain (constant x) ---

# 1. Find the middle x value (center of the channel)
x_unique = np.unique(np.round(x, 3))
x_center = x_unique[len(x_unique)//2]
tol = (x_unique[1] - x_unique[0]) / 2 if len(x_unique) > 1 else 1.0

# 2. Select indices where x is close to x_center
cross_idx = np.where(np.abs(x - x_center) < tol)[0]
y_cross = y[cross_idx]
wd_cross = wd[t_index, cross_idx]

# 3. Analytical solution for waterdepth(y)
T_sidereal = 23*3600 + 56*60 + 4.1  # seconds
Omega = 2 * np.pi / T_sidereal
phi = np.deg2rad(45)
f_coriolis = 2 * Omega * np.sin(phi)
y0 = 150000.0
g = 9.81

if u is not None:
    u_cross = u[t_index, cross_idx]
    u_mean = np.mean(u_cross)
else:
    u_mean = 1.0  # fallback value

wd_analytical_cross = -f_coriolis * (y_cross - y0) / g * u_mean

# Get bed elevation if available
if 'mesh2d_face_z' in data_vars:
    bed_z = data_vars['mesh2d_face_z'][:]
    bed_cross = bed_z[cross_idx]
else:
    bed_cross = np.zeros_like(wd_cross)  # fallback if bed elevation is not present

# Calculate water elevation
elev_cross = wd_cross + bed_cross

# Analytical water elevation (assuming analytical water depth is relative to same bed)
elev_analytical_cross = wd_analytical_cross + bed_cross

# 4. Sort by y for plotting
sort_y_idx = np.argsort(y_cross)
y_cross_sorted = y_cross[sort_y_idx]
elev_cross_sorted = elev_cross[sort_y_idx]
elev_analytical_sorted = elev_analytical_cross[sort_y_idx]

# 5. Plot
plt.figure(figsize=(9, 5))
plt.plot(y_cross_sorted, elev_cross_sorted-500, 'bo-', label='Numerical water elevation')
plt.plot(y_cross_sorted, elev_analytical_sorted, 'r--', label='Analytical water elevation (Coriolis)')
plt.xlabel('y [m]')
plt.ylabel('Water elevation [m]')
plt.title(f'Cross-section at x = {x_center:.1f} m, t = {t:.2f} s')
plt.legend()
plt.grid(True, linestyle=':')
plt.tight_layout()
plt.show()
