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
MDU_FILE_NAME = "planar1d"
MAP_FILE_NAME = "planar1d_map.nc"
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
wd = data_vars['mesh2d_waterdepth'][:]  # water depth (time, face)
vel = data_vars['mesh2d_ucmag'][:] if 'mesh2d_ucmag' in data_vars else None  # velocity magnitude (time, face)
times = data_vars['time'][:] if 'time' in data_vars else np.arange(wd.shape[0])

# --- Analytical solution functions ---
def analytical_waterlevel(x, t, eta0, h0, r0, omega):
    # h(x, t) = (eta0 * h0 / r0) * (2x cos(omega t) - eta0 r0 cos^2(omega t))
    return (eta0 * h0 / r0) * (2 * (x) * np.cos(omega * t) - eta0 * r0 * np.cos(omega * t)**2)

def analytical_velocity(x, t, eta0, r0, omega):
    # u(x, t) = -eta0 * r0 * omega * sin(omega * t)
    return abs(-eta0 * r0 * omega * np.sin(omega * t) * np.ones_like(x))

# --- Parameters for analytical solution ---
eta0 = 0.23   # Amplitude [m]
h0 = 10.0  # Reference water depth [m]
r0 = 120.0   # Reference length [m]
omega = 2*np.pi/53.8284  # Angular frequency [1/s], e.g. period = 54s

# Add this setting at the top of your plotting section:
PLOT_COMBINED = True  # Set to False to plot individually, True for combined plot
PLOT_TEMPORAL_EVOLUTION = True  # Set to True to enable the temporal evolution plot

if PLOT_COMBINED:
    # --- Combined Water Depth, Velocity, and Water Level/Bed Level Plot ---
    plt.close('all')  # Close all existing plots before creating a new one

    # Prepare data
    x = data_vars['mesh1d_node_x'][:]  # x-coordinates of nodes
    wd = data_vars['mesh1d_waterdepth'][:]  # water depth (time, node)
    vel = data_vars['mesh1d_ucmag'][:] if 'mesh1d_ucmag' in data_vars else None  # velocity magnitude (time, node)
    wl = data_vars['mesh1d_s1'][:] if 'mesh1d_s1' in data_vars else None  # water level (time, node)
    bed = wl - wd if wl is not None else None  # bed level (time, node)
    times = data_vars['time'][:] if 'time' in data_vars else np.arange(wd.shape[0])

    # Create figure with three subplots side by side
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 6), sharex=True)
    plt.subplots_adjust(bottom=0.22, wspace=0.35)  # More space for slider and between plots

    # Water depth plot (left)
    y_min_wd = 0.9 * np.nanmin(wd)
    y_max_wd = 1.1 * np.nanmax(wd)
    ax1.set_ylim(y_min_wd, y_max_wd)
    line1, = ax1.plot(x, wd[0, :] if wd.ndim == 2 else wd[:], 'b-', linewidth=lw, label='Water depth')
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('Water depth [m]')
    ax1.set_title(f'Water depth along channel (t = {times[0]:.2f})')
    ax1.grid(True, linestyle=':')
    ax1.legend()

    # Velocity plot (middle)
    if vel is not None:
        y_min_vel = 0.9 * np.nanmin(vel)
        y_max_vel = 1.1 * np.nanmax(vel)
        ax2.set_ylim(y_min_vel, y_max_vel)
        line2, = ax2.plot(x, vel[0, :] if vel.ndim == 2 else vel[:], 'r-', linewidth=lw, label='Velocity magnitude')
        ax2.set_xlabel('x [m]')
        ax2.set_ylabel('Velocity [m/s]')
        ax2.set_title(f'Velocity magnitude along channel (t = {times[0]:.2f})')
        ax2.grid(True, linestyle=':')
        ax2.legend()
    else:
        line2 = None
        ax2.set_visible(False)
        print("Variable 'mesh2d_ucmag' not found in NetCDF file for velocity plot.")

    # Water level and bed level plot (right)
    if wl is not None and bed is not None:
        y_min_wl = min(np.nanmin(wl), np.nanmin(bed))
        y_max_wl = max(np.nanmax(wl), np.nanmax(bed))
        ax3.set_ylim(y_min_wl-0.1*y_max_wl, 1.1 * y_max_wl)
        line_wl, = ax3.plot(x, wl[0, :] if wl.ndim == 2 else wl[:], 'b-', linewidth=lw, label='Water level')
        line_bed, = ax3.plot(x, bed[0, :] if bed.ndim == 2 else bed[:], 'k-', linewidth=lw, label='Bed level')
        ax3.set_xlabel('x [m]')
        ax3.set_ylabel('Level [m]')
        ax3.set_title(f'Water level and bed level (t = {times[0]:.2f})')
        ax3.grid(True, linestyle=':')
        ax3.legend()
    else:
        line_wl = line_bed = None
        ax3.set_visible(False)
        print("Variable 'mesh2d_s1' (water level) not found in NetCDF file")

    plt.tight_layout(rect=[0, 0.12, 1, 1])  # Leave space for slider

    # Slider axis (universal for all subplots)
    ax_slider = plt.axes([0.15, 0.07, 0.7, 0.03])  # [left, bottom, width, height]
    slider = Slider(ax_slider, 'Time step', 0, wd.shape[0]-1, valinit=0, valstep=1)

    def update(val):
        idx = int(slider.val)
        t = times[idx]
        line1.set_ydata(wd[idx, :] if wd.ndim == 2 else wd[:])
        ax1.set_title(f'Water depth along channel (t = {t:.2f} s)')
        if line2 is not None:
            line2.set_ydata(vel[idx, :] if vel.ndim == 2 else vel[:])
            ax2.set_title(f'Velocity magnitude along channel (t = {t:.2f} s)')
        if line_wl is not None and line_bed is not None:
            line_wl.set_ydata(wl[idx, :] if wl.ndim == 2 else wl[:])
            line_bed.set_ydata(bed[idx, :] if bed.ndim == 2 else bed[:])
            ax3.set_title(f'Water level and bed level (t = {t:.2f} s)')

        # --- Analytical overlays ---
        ana_wl = analytical_waterlevel(x, t, eta0, h0, r0, omega)
        ana_vel = analytical_velocity(x, t, eta0, r0, omega)
        # Plot or update analytical overlays
        if not hasattr(update, "ana_wl_line"):
            update.ana_wl_line, = ax3.plot(x+150, ana_wl, 'g--', label='Analytical WL')
            ax3.legend()
        else:
            update.ana_wl_line.set_ydata(ana_wl)
        if line2 is not None:
            if not hasattr(update, "ana_vel_line"):
                update.ana_vel_line, = ax2.plot(x, ana_vel, 'g--', label='Analytical vel')
                ax2.legend()
            else:
                update.ana_vel_line.set_ydata(ana_vel)

        fig.canvas.draw_idle()

    # Initial analytical overlays
    t0 = times[0]
    ana_wl0 = analytical_waterlevel(x, t0, eta0, h0, r0, omega)
    ana_vel0 = analytical_velocity(x, t0, eta0, r0, omega)
    ana_wd_line, = ax1.plot(x, ana_wl0 - bed[0, :]-5.78, 'g--', label='Analytical WD')
    if line2 is not None:
        ana_vel_line, = ax2.plot(x, ana_vel0, 'g--', label='Analytical vel')
    ax1.legend()
    if line2 is not None:
        ax2.legend()

    slider.on_changed(update)
    plt.show()

else:
    # --- Water depth plot (individual) ---
    plt.close('all')
    fig, ax = plt.subplots(figsize=(10, 5))
    plt.subplots_adjust(bottom=0.22)
    y_min = 0.9 * np.nanmin(wd)
    y_max = 1.1 * np.nanmax(wd)
    ax.set_ylim(y_min, y_max)
    line, = ax.plot(x, wd[0, :] if wd.ndim == 2 else wd[:], 'b-', linewidth=lw, label='Water depth')
    ax.set_xlabel('x [m]', fontsize=fs)
    ax.set_ylabel('Water depth [m]', fontsize=fs)
    ax.set_title(f'Water depth along channel (t = {times[0]:.2f})', fontsize=fs)
    ax.grid(True, linestyle=':')
    ax.legend()
    plt.tight_layout()
    ax_slider = plt.axes([0.15, 0.07, 0.7, 0.03])
    slider = Slider(ax_slider, 'Time step', 0, wd.shape[0]-1, valinit=0, valstep=1)

    def update(val):
        idx = int(slider.val)
        line.set_ydata(wd[idx, :] if wd.ndim == 2 else wd[:])
        ax.set_title(f'Water depth along channel (t = {times[idx]:.2f})', fontsize=fs)
        fig.canvas.draw_idle()

    slider.on_changed(update)
    plt.show()

    # Save water depth figure at final timestep
    final_idx = wd.shape[0] - 1
    line.set_ydata(wd[final_idx, :] if wd.ndim == 2 else wd[:])
    ax.set_title(f'Water depth along channel (t = {times[final_idx]:.2f})', fontsize=fs)
    fig.savefig(os.path.join(FIGURES, "water_depth_final.png"))

    # --- Velocity plot (individual) ---
    if vel is not None:
        fig2, ax2 = plt.subplots(figsize=(10, 5))
        plt.subplots_adjust(bottom=0.22)
        y_min_v = 0.9 * np.nanmin(vel)
        y_max_v = 1.1 * np.nanmax(vel)
        ax2.set_ylim(y_min_v, y_max_v)
        line2, = ax2.plot(x, vel[0, :] if vel.ndim == 2 else vel[:], 'r-', linewidth=lw, label='Velocity magnitude')
        ax2.set_xlabel('x [m]', fontsize=fs)
        ax2.set_ylabel('Velocity [m/s]', fontsize=fs)
        ax2.set_title(f'Velocity magnitude along channel (t = {times[0]:.2f})', fontsize=fs)
        ax2.grid(True, linestyle=':')
        ax2.legend()
        plt.tight_layout()
        ax_slider2 = plt.axes([0.15, 0.07, 0.7, 0.03])
        slider2 = Slider(ax_slider2, 'Time step', 0, vel.shape[0]-1, valinit=0, valstep=1)

        def update2(val):
            idx = int(slider2.val)
            line2.set_ydata(vel[idx, :] if vel.ndim == 2 else vel[:])
            ax2.set_title(f'Velocity magnitude along channel (t = {times[idx]:.2f})', fontsize=fs)
            fig2.canvas.draw_idle()

        slider2.on_changed(update2)
        plt.show()

        # Save velocity figure at final timestep
        final_idx2 = vel.shape[0] - 1
        line2.set_ydata(vel[final_idx2, :] if vel.ndim == 2 else vel[:])
        ax2.set_title(f'Velocity magnitude along channel (t = {times[final_idx2]:.2f})', fontsize=fs)
        fig2.savefig(os.path.join(FIGURES, "velocity_final.png"))
    else:
        print("Variable 'mesh2d_ucmag' not found in NetCDF file")

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

    # Use mesh1d_node_x if available, otherwise mesh2d_face_x
    if 'mesh1d_node_x' in data_vars:
        x_arr = data_vars['mesh1d_node_x'][:]
        wl_arr = data_vars['mesh1d_s1'][:] if 'mesh1d_s1' in data_vars else None
        vel_arr = data_vars['mesh1d_ucmag'][:] if 'mesh1d_ucmag' in data_vars else None
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
    ana_wl_series = np.array([analytical_waterlevel(x_actual-150, t, eta0, h0, r0, omega) for t in times_plot])
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