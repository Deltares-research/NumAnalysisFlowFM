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
MDU_FILE_NAME = "radial2d"
MAP_FILE_NAME = "radial2d_map.nc"
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
def analytical_waterlevel(r, t):
    denom = 1 - A * np.cos(omega * t)
    sqrt_term = np.sqrt(1 - A**2)
    term1 = D0 * (sqrt_term / denom - 1)
    term2 = (r**2 / L**2) * ((1 - A**2) / denom**2 - 1)
    return term1 - D0 * term2

def analytical_velocity(r, t):
    denom = 1 - A * np.cos(omega * t)
    return 0.5 * r * omega * (A * np.sin(omega * t)) / denom


# Add this setting at the top of your plotting section:
PLOT_COMBINED = True  # Set to False to plot individually, True for combined plot
PLOT_TEMPORAL_EVOLUTION = False  # Set to True to enable the temporal evolution plot

from matplotlib.widgets import Slider

if PLOT_COMBINED:
    plt.close('all')

    # --- Load mesh and data --- 
    x = np.asarray(data_vars['mesh2d_face_x'][:])
    y = np.asarray(data_vars['mesh2d_face_y'][:])
    r = np.sqrt(x**2 + y**2)

    # Handle masked arrays: fill masked values with 0.0
    wd_raw = data_vars['mesh2d_waterdepth'][:]
    wd = np.ma.filled(wd_raw, 0.0)

    wl_raw = data_vars.get('mesh2d_s1', None)
    wl = np.ma.filled(wl_raw, 0.0) if wl_raw is not None else None

    u_raw = data_vars.get('mesh2d_ucx', None)
    v_raw = data_vars.get('mesh2d_ucy', None)
    u = np.ma.filled(u_raw, 0.0) if u_raw is not None else None
    v = np.ma.filled(v_raw, 0.0) if v_raw is not None else None

    # Compute bed elevation = water level - depth (t=1), safely
    bed = wl[1] - wd[1] if wl is not None else None

    # Time axis
    times = data_vars['time'][:] if 'time' in data_vars else np.arange(wd.shape[0])
    n_times = wd.shape[0]


    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 6), sharex=True)
    plt.subplots_adjust(bottom=0.22, wspace=0.35)

    sort_idx = np.argsort(r)
    r_sorted = r[sort_idx]

    # --- Initial plots ---
    tri_contour1 = ax1.tricontourf(x, y, wd[0]-bed[:], levels=20, cmap='viridis')
    cbar1 = fig.colorbar(tri_contour1, ax=ax1, label='Water level [m]')
    ax1.set_title(f'Water level (t = {times[0]:.2f} s)')
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')
    ax1.grid(True, linestyle=':')

    if u is not None and v is not None:
        vel = np.sqrt(u**2 + v**2)
        tri_contour2 = ax2.tricontourf(x,y, vel[0], levels=20, cmap='plasma')
        cbar2 = fig.colorbar(tri_contour2, ax=ax2, label='Velocity [m/s]')
        ax2.set_title(f'Velocity magnitude (t = {times[0]:.2f} s)')
        ax2.set_xlabel('r [m]')
        ax2.set_ylabel('Velocity [m/s]')
        ax2.grid(True, linestyle=':')
    else:
        vel_line = None
        print("Velocity data not found.")

    ax3.axis('off')  # placeholder

    ax_slider = plt.axes([0.25, 0.08, 0.5, 0.03])
    slider = Slider(ax_slider, 'Timestep', 0, n_times - 1, valinit=0, valstep=1)

    # --- Store colorbars so we can remove them ---
    colorbars = [cbar1, cbar2 if u is not None and v is not None else None]

    def update(val):
        idx = int(slider.val)

        # Remove previous contours and colorbars
        if colorbars[0] is not None:
            colorbars[0].remove()

        tri_contour1 = ax1.tricontourf(x, y, wd[idx]-bed[:], levels=20, cmap='viridis')
        colorbars[0] = fig.colorbar(tri_contour1, ax=ax1, label='Water level [m]')
        ax1.set_title(f'Water depth (t = {times[idx]:.2f} s)')

        if u is not None and v is not None:
            if colorbars[1] is not None:
                colorbars[1].remove()
            tri_contour2 = ax2.tricontourf(x, y, vel[idx], levels=20, cmap='plasma')
            colorbars[1] = fig.colorbar(tri_contour2, ax=ax2, label='Velocity [m/s]')
            ax2.set_title(f'Velocity magnitude (t = {times[idx]:.2f} s)')

        fig.canvas.draw_idle()

    slider.on_changed(update)
    plt.show()


    # Calculate water level and bed level
    if 'mesh2d_s1' in data_vars:
        wl = data_vars['mesh2d_s1'][:]  # water level (time, face)
        bed = wl - wd                   # bed level (time, face)
    
        # Plot water level and bed level in the same figure
        plt.close('all')
        fig3, ax3 = plt.subplots(figsize=(10, 5))
        plt.subplots_adjust(bottom=0.22)
        # Plot initial timestep
        line_wl, = ax3.plot(x, wl[0, :] if wl.ndim == 2 else wl[:], 'b-', linewidth=lw, label='Water level')
        line_bed, = ax3.plot(x, bed[0, :] if bed.ndim == 2 else bed[:], 'k-', linewidth=lw, label='Bed level')
        ax3.set_xlabel('x [m]', fontsize=fs)
        ax3.set_ylabel('Level [m]', fontsize=fs)
        ax3.set_title(f'Water level and bed level along channel (t = {times[0]:.2f})', fontsize=fs)
        ax3.grid(True, linestyle=':')
        ax3.legend()
        plt.tight_layout()
        # Slider for this plot
        ax_slider3 = plt.axes([0.15, 0.07, 0.7, 0.03])
        slider3 = Slider(ax_slider3, 'Time step', 0, wl.shape[0]-1, valinit=0, valstep=1)
    
        def update3(val):
            idx = int(slider3.val)
            line_wl.set_ydata(wl[idx, :] if wl.ndim == 2 else wl[:])
            line_bed.set_ydata(bed[idx, :] if bed.ndim == 2 else bed[:])
            ax3.set_title(f'Water level and bed level along channel (t = {times[idx]:.2f})', fontsize=fs)
            fig3.canvas.draw_idle()
    
        slider3.on_changed(update3)
        plt.show()
    
        # Save figure at final timestep
        final_idx3 = wl.shape[0] - 1
        line_wl.set_ydata(wl[final_idx3, :] if wl.ndim == 2 else wl[:])
        line_bed.set_ydata(bed[final_idx3, :] if bed.ndim == 2 else bed[:])
        ax3.set_title(f'Water level and bed level along channel (t = {times[final_idx3]:.2f})', fontsize=fs)
        fig3.savefig(os.path.join(FIGURES, "waterlevel_bedlevel_final.png"))
    else:
        print("Variable 'mesh2d_s1' (water level) not found in NetCDF file")

if PLOT_TEMPORAL_EVOLUTION:
    # --- Temporal evolution at a specific x-location ---
    # Customize the x-location here (in meters)
    x_probe = 200  # Change this value as needed

    # Use mesh2d_node_x if available, otherwise mesh2d_face_x
    if 'mesh2d_node_x' in data_vars:
        x_arr = data_vars['mesh2d_node_x'][:]
        wl_arr = data_vars['mesh2d_s1'][:] if 'mesh2d_s1' in data_vars else None
        vel_arr = data_vars['mesh2d_ucmag'][:] if 'mesh2d_ucmag' in data_vars else None
    elif 'mesh2d_face_x' in data_vars:
        x_arr = data_vars['mesh2d_face_x'][:]
        wl_arr = data_vars['mesh2d_s1'][:] if 'mesh2d_s1' in data_vars else None
        vel_arr = data_vars['mesh2d_ucmag'][:] if 'mesh2d_ucmag' in data_vars else None
    else:
        x_arr = x
        wl_arr = wl
        vel_arr = vel

    # Find the index of the node/face closest to x_probe
    idx_probe = np.argmin(np.abs(x_arr - x_probe))
    x_actual = x_arr[idx_probe]
    print(f"Temporal evolution will be shown at x = {x_actual:.2f} m (closest to requested {x_probe} m)")

    # Prepare time series
    times_plot = times
    wl_series = wl_arr[:, idx_probe] if wl_arr is not None else None
    vel_series = vel_arr[:, idx_probe] if vel_arr is not None else None

    # Analytical solutions at this location
    ana_wl_series = np.array([analytical_waterlevel(x_actual t, eta0, h0, r0, omega) for t in times_plot])
    ana_vel_series = np.array([analytical_velocity(x_actual, t, eta0, r0, omega) for t in times_plot])

    # Plot
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
    # Water level (left)
    if wl_series is not None:
        axL.plot(times_plot, wl_series, 'b-', label='Numerical WL')
    axL.plot(times_plot, ana_wl_series, 'g--', label='Analytical WL')
    axL.set_xlabel('Time [s]')
    axL.set_ylabel('Water level [m]')
    axL.set_title(f'Water level at x = {x_actual:.2f} m')
    axL.legend()
    axL.grid(True, linestyle=':')

    # Velocity (right)
    if vel_series is not None:
        axR.plot(times_plot, vel_series, 'r-', label='Numerical vel')
    axR.plot(times_plot, ana_vel_series, 'g--', label='Analytical vel')
    axR.set_xlabel('Time [s]')
    axR.set_ylabel('Velocity [m/s]')
    axR.set_title(f'Velocity at x = {x_actual:.2f} m')
    axR.legend()
    axR.grid(True, linestyle=':')

    plt.tight_layout()
    plt.show()


    import matplotlib.pyplot as plt


# Time setup
times = data_vars['time'][:] if 'time' in data_vars else np.arange(wd.shape[0])
t_index = 5
t = times[t_index]



# --- Evaluate analytical solution ---
wl_analytical = analytical_waterlevel(r, t)

# --- Handle mask in numerical water level ---
wl_num = np.ma.filled(wl[t_index], 0.0)

# --- Compute difference ---
diff = wl_num - wl_analytical

# --- Optional: mask invalid points (e.g. outside domain)
if np.ma.is_masked(wl[t_index]):
    diff = np.ma.array(diff, mask=wl[t_index].mask)
    
# --- Compute relative difference ---
# Avoid division by zero or near-zero values
tol = 1e-6
mask = np.abs(wl_analytical) < tol
relative_diff = np.zeros_like(wl_analytical)
relative_diff[~mask] = (wl_num[~mask] - wl_analytical[~mask]) / wl_analytical[~mask]
relative_diff[mask] = 0.0  # or np.nan, depending on what you prefer

# --- Optional: apply original mask as well
if np.ma.is_masked(wl[t_index]):
    relative_diff = np.ma.array(relative_diff, mask=wl[t_index].mask | mask)

# --- Plotting ---
plt.figure(figsize=(10, 5))
sc = plt.tripcolor(x, y, relative_diff, shading='flat', cmap='RdBu',
                   vmin=-np.max(np.abs(relative_diff)), vmax=np.max(np.abs(relative_diff)))
plt.colorbar(sc, label="Relative Water Level Error\n(numerical - analytical) / analytical")
plt.title(f"Relative Water Level Error at t = {t:.2f} s (Index {t_index})")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.axis("equal")
plt.tight_layout()
plt.show()
