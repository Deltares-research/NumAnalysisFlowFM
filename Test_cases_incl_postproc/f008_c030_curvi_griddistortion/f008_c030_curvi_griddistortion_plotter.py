"""This script is a general python script to plot results of 1D cases that are performed in 2D using D-Hydro FM"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import os
import configparser
from matplotlib.ticker import MaxNLocator, MultipleLocator
from matplotlib.widgets import Slider
import matplotlib as mpl
try:
    from scipy.interpolate import griddata
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("Warning: scipy not available. Using alternative interpolation method.")

# Set global font to Times New Roman, size 12, and mathtext italic
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.size'] = 12
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Times New Roman'
mpl.rcParams['mathtext.it'] = 'Times New Roman:italic'
mpl.rcParams['axes.titlesize'] = 12
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['legend.fontsize'] = 12
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12

# Set working directory

PATH = os.getcwd() # Ensure that this file is present in the test folder that contains ...\dflowfmoutput
OUTPUT_FOLDER = "dflowfmoutput" # The folder where the output files are stored
MDU_FILE_NAME = "curvidistortion"
MAP_FILE_NAME = "curvidistortion_map.nc"
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
s1 = data_vars['mesh2d_s1'][:]  # water level (time, face) 
bl = data_vars['mesh2d_flowelem_bl'][:]  # bed level (face)
vel = data_vars['mesh2d_ucmag'][:] if 'mesh2d_ucmag' in data_vars else None  # velocity magnitude (time, face)
times = data_vars['time'][:] if 'time' in data_vars else np.arange(wd.shape[0])

# Use the last time step for plotting (index -1)
time_idx = -1

# Debug: Let's examine the water depth data more closely
print("\n=== DATA ANALYSIS ===")
print(f"Water depth array shape: {wd.shape}")
print(f"Water depth (wd) range: {np.min(wd):.3f} to {np.max(wd):.3f} m")
print(f"Water level (s1) range: {np.min(s1):.3f} to {np.max(s1):.3f} m")
print(f"Bed level (bl) range: {np.min(bl):.3f} to {np.max(bl):.3f} m")
print(f"Calculated depth (s1-bl) range: {np.min(s1-bl):.3f} to {np.max(s1-bl):.3f} m")

# Check if water depth matches s1-bl
depth_diff = np.abs(wd - (s1 - bl))
print(f"Max difference between wd and (s1-bl): {np.max(depth_diff):.6f} m")

# Velocity analysis
if vel is not None:
    print(f"\nVelocity magnitude (ucmag) analysis:")
    print(f"  Shape: {vel.shape}")
    print(f"  Range: {np.min(vel):.6f} to {np.max(vel):.6f} m/s")
    print(f"  Mean: {np.mean(vel):.6f} m/s")
    print(f"  Std: {np.std(vel):.6f} m/s")
    print(f"  Non-zero values: {np.sum(vel > 1e-10)}")
    print(f"  Zero/near-zero values: {np.sum(vel <= 1e-10)}")
    
    # Check individual velocity components
    ucx = data_vars['mesh2d_ucx'][:]  # x-component
    ucy = data_vars['mesh2d_ucy'][:]  # y-component
    print(f"\nVelocity components:")
    print(f"  ucx range: {np.min(ucx):.6f} to {np.max(ucx):.6f} m/s")
    print(f"  ucy range: {np.min(ucy):.6f} to {np.max(ucy):.6f} m/s")
    
    # Check if velocity magnitude matches sqrt(ucx^2 + ucy^2)
    vel_calculated = np.sqrt(ucx**2 + ucy**2)
    vel_diff = np.abs(vel - vel_calculated)
    print(f"  Max difference between ucmag and sqrt(ucx^2+ucy^2): {np.max(vel_diff):.6e}")
    
    # For final time step
    print(f"\nVelocity at time step {time_idx}:")
    print(f"  ucmag range: {np.min(vel[time_idx, :]):.6f} to {np.max(vel[time_idx, :]):.6f} m/s")
    print(f"  ucx range: {np.min(ucx[time_idx, :]):.6f} to {np.max(ucx[time_idx, :]):.6f} m/s")
    print(f"  ucy range: {np.min(ucy[time_idx, :]):.6f} to {np.max(ucy[time_idx, :]):.6f} m/s")
    
    # Check for unrealistic values
    high_vel_mask = vel[time_idx, :] > 10.0  # velocities > 10 m/s
    if np.any(high_vel_mask):
        print(f"  WARNING: {np.sum(high_vel_mask)} faces have velocities > 10 m/s")
        print(f"  Max velocity location indices: {np.where(high_vel_mask)[0][:5]}")  # Show first 5

# Let's check what the expected velocity should be for this type of flow
print(f"Expected velocity for steady shallow water flow should be << 10 m/s")

# Let's check what the actual velocity data contains

# Let's check what the actual analytical depth should be
print(f"\nActual water depth values for time step {time_idx}:")
print(f"  Mean: {np.mean(wd[time_idx, :]):.3f} m")
print(f"  Std:  {np.std(wd[time_idx, :]):.3f} m")
print(f"  Min:  {np.min(wd[time_idx, :]):.3f} m")
print(f"  Max:  {np.max(wd[time_idx, :]):.3f} m")

# Check if the data is actually varying spatially or if there's something else
print(f"\nWater level (s1) for time step {time_idx}:")
print(f"  Mean: {np.mean(s1[time_idx, :]):.3f} m")
print(f"  Std:  {np.std(s1[time_idx, :]):.3f} m")
print(f"  Min:  {np.min(s1[time_idx, :]):.3f} m")
print(f"  Max:  {np.max(s1[time_idx, :]):.3f} m")

print(f"\nBed level (bl):")
print(f"  Mean: {np.mean(bl):.3f} m")
print(f"  Std:  {np.std(bl):.3f} m")
print(f"  Min:  {np.min(bl):.3f} m")
print(f"  Max:  {np.max(bl):.3f} m")
print("===================\n")

# Split the dataset into 2 grids based on y-coordinates
y_min = np.min(y)
y_max = np.max(y)
y_mid = y_max / 2

# Create indices for each grid
grid_indices = []
# Grid 1: fully below y_max/2
mask1 = y < y_mid
grid_indices.append(mask1)

# Grid 2: fully above y_max/2
mask2 = y >= y_mid
grid_indices.append(mask2)

print(f"Y range: {y_min:.3f} to {y_max:.3f}")
print(f"Y midpoint: {y_mid:.3f}")
for i, mask in enumerate(grid_indices):
    print(f"Grid {i+1}: {np.sum(mask)} faces")

# Create figure with 2 rows and 2 columns
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Water Depth (Constant) and Velocity (Constant) Distribution Across 2 Grid Levels', fontsize=14, fontweight='bold')

for grid_num in range(2):
    mask = grid_indices[grid_num]
    
    if np.sum(mask) == 0:
        print(f"Warning: Grid {grid_num+1} has no data points")
        continue
    
    # Extract data for this grid
    x_grid = x[mask]
    y_grid = y[mask]
    wd_grid = wd[time_idx, mask]
    
    # Water depth plot (left column) - Since water depth is constant, use scatter plot
    ax_wd = axes[grid_num, 0]
    
    # Since water depth is constant, scatter plot is more appropriate than contour
    # Use actual data range for colorbar
    wd_min, wd_max = np.min(wd_grid), np.max(wd_grid)
    scatter_wd = ax_wd.scatter(x_grid, y_grid, c=wd_grid, cmap='Blues', s=50, alpha=0.8, 
                              vmin=wd_min, vmax=wd_max)
    
    ax_wd.set_title(f'Grid {grid_num+1}: Water Depth (Constant)', fontweight='bold')
    ax_wd.set_xlabel('X Coordinate (m)')
    ax_wd.set_ylabel('Y Coordinate (m)')
    ax_wd.grid(True, alpha=0.3)
    
    # Add colorbar for water depth
    cbar_wd = plt.colorbar(scatter_wd, ax=ax_wd, shrink=0.8)
    cbar_wd.set_label('Water Depth (m)', rotation=270, labelpad=20)
    
    # Velocity plot (right column)
    ax_vel = axes[grid_num, 1]
    if vel is not None:
        vel_grid = vel[time_idx, mask]
        
        # Check if velocity is approximately constant (like water depth)
        vel_std = np.std(vel_grid)
        vel_mean = np.mean(vel_grid)
        
        if vel_std < 0.01 * vel_mean:  # If velocity is nearly constant (1% variation)
            # Use scatter plot for constant/nearly constant velocity
            # Use actual data range for colorbar
            vel_min, vel_max = np.min(vel_grid), np.max(vel_grid)
            scatter_vel = ax_vel.scatter(x_grid, y_grid, c=vel_grid, cmap='Reds', s=50, alpha=0.8,
                                       vmin=vel_min, vmax=vel_max)
            contour_vel = scatter_vel
            print(f"Grid {grid_num+1}: Velocity is nearly constant ({vel_mean:.3f} ± {vel_std:.6f} m/s), using scatter plot")
        else:
            # Check if we have enough points for triangulation
            if len(x_grid) >= 3:
                try:
                    if HAS_SCIPY:
                        # Create regular grid for interpolation for velocity
                        x_min, x_max = x_grid.min(), x_grid.max()
                        y_min_grid, y_max_grid = y_grid.min(), y_grid.max()
                        
                        # Create grid points for interpolation
                        grid_resolution = 50  # Adjust for finer/coarser resolution
                        xi = np.linspace(x_min, x_max, grid_resolution)
                        yi = np.linspace(y_min_grid, y_max_grid, grid_resolution)
                        xi_mesh, yi_mesh = np.meshgrid(xi, yi)
                        
                        # Interpolate velocity data onto regular grid
                        vel_interp = griddata((x_grid, y_grid), vel_grid, (xi_mesh, yi_mesh), method='linear')
                        
                        # Create filled contour plot for velocity with actual data range
                        vel_min, vel_max = np.min(vel_grid), np.max(vel_grid)
                        contour_vel = ax_vel.contourf(xi_mesh, yi_mesh, vel_interp, levels=20, cmap='Reds', alpha=0.8,
                                                     vmin=vel_min, vmax=vel_max)
                    else:
                        # Use triangular contour plot with actual data range
                        vel_min, vel_max = np.min(vel_grid), np.max(vel_grid)
                        contour_vel = ax_vel.tricontourf(x_grid, y_grid, vel_grid, levels=20, cmap='Reds', alpha=0.8,
                                                        vmin=vel_min, vmax=vel_max)
                except (RuntimeError, ValueError) as e:
                    # Fallback to scatter plot if triangulation fails
                    print(f"Triangulation failed for velocity in grid {grid_num+1}, using scatter plot: {e}")
                    scatter_vel = ax_vel.scatter(x_grid, y_grid, c=vel_grid, cmap='Reds', s=50, alpha=0.8)
                    contour_vel = scatter_vel
            else:
                # Not enough points for contour, use scatter plot
                print(f"Not enough points for velocity contour in grid {grid_num+1}, using scatter plot")
                scatter_vel = ax_vel.scatter(x_grid, y_grid, c=vel_grid, cmap='Reds', s=50, alpha=0.8)
                contour_vel = scatter_vel
        
        ax_vel.set_title(f'Grid {grid_num+1}: Velocity Magnitude', fontweight='bold')
        
        # Add colorbar for velocity
        cbar_vel = plt.colorbar(contour_vel, ax=ax_vel, shrink=0.8)
        cbar_vel.set_label('Velocity Magnitude (m/s)', rotation=270, labelpad=20)
    else:
        ax_vel.text(0.5, 0.5, 'Velocity data\nnot available', 
                   ha='center', va='center', transform=ax_vel.transAxes,
                   fontsize=12, bbox=dict(boxstyle='round', facecolor='lightgray'))
        ax_vel.set_title(f'Grid {grid_num+1}: Velocity Magnitude', fontweight='bold')
    
    ax_vel.set_xlabel('X Coordinate (m)')
    ax_vel.set_ylabel('Y Coordinate (m)')
    ax_vel.grid(True, alpha=0.3)

# Adjust layout to prevent overlap
plt.tight_layout()

# Create Figures directory if it doesn't exist
if not os.path.exists(FIGURES):
    os.makedirs(FIGURES)

# Save the figure
plt.savefig(os.path.join(FIGURES, 'grid_analysis_contours.png'), dpi=300, bbox_inches='tight')
plt.show()

# Create a second figure for water depth comparison with analytical solution
analytical_depth = 2.037  # Constant analytical water depth in meters

fig2, ax_comparison = plt.subplots(1, 1, figsize=(12, 8))

# Plot water depth for each grid with different colors
colors = ['blue', 'red']
grid_labels = ['Grid 1', 'Grid 2']

for grid_num in range(2):
    mask = grid_indices[grid_num]
    
    if np.sum(mask) == 0:
        continue
    
    # Extract data for this grid
    x_grid = x[mask]
    wd_grid = wd[time_idx, mask]
    
    # Debug: print stats for each grid
    print(f"Grid {grid_num+1}: wd range {np.min(wd_grid):.3f} to {np.max(wd_grid):.3f} m")
    
    # Plot water depth vs x-coordinate
    ax_comparison.scatter(x_grid, wd_grid, c=colors[grid_num], label=grid_labels[grid_num], 
                         alpha=0.7, s=20)

# Plot analytical solution as a horizontal line
x_range = [np.min(x), np.max(x)]
ax_comparison.plot(x_range, [analytical_depth, analytical_depth], 'k--', linewidth=2, 
                  label=f'Analytical Solution (h = {analytical_depth} m)')

ax_comparison.set_xlabel('X Coordinate (m)')
ax_comparison.set_ylabel('Water Depth (m)')
ax_comparison.set_title('Water Depth Comparison: Numerical vs Analytical Solution', fontweight='bold')
ax_comparison.grid(True, alpha=0.3)
ax_comparison.legend()

# Set y-axis limits to better show the comparison
y_min_plot = min(np.min(wd[time_idx, :]), analytical_depth) * 0.95
y_max_plot = max(np.max(wd[time_idx, :]), analytical_depth) * 1.05
ax_comparison.set_ylim(y_min_plot, y_max_plot)

plt.tight_layout()

# Save the comparison figure
plt.savefig(os.path.join(FIGURES, 'water_depth_comparison.png'), dpi=300, bbox_inches='tight')
plt.show()

print(f"Contour plots saved to: {os.path.join(FIGURES, 'grid_analysis_contours.png')}")
print(f"Comparison plot saved to: {os.path.join(FIGURES, 'water_depth_comparison.png')}")
print(f"Analytical water depth: {analytical_depth} m")
print(f"Numerical water depth range: {np.min(wd[time_idx, :]):.3f} to {np.max(wd[time_idx, :]):.3f} m")
print(f"Difference from analytical: {abs(np.mean(wd[time_idx, :]) - analytical_depth):.6f} m")

# Create a third figure for convergence analysis
fig3, ax_convergence = plt.subplots(1, 1, figsize=(10, 8))

# Calculate L2 norm (RMS) for each grid level
grid_levels = []
l2_norms = []
grid_sizes = []

for grid_num in range(2):
    mask = grid_indices[grid_num]
    
    if np.sum(mask) == 0:
        continue
    
    # Extract data for this grid
    wd_grid = wd[time_idx, mask]
    n_faces = len(wd_grid)
    
    # Calculate differences from analytical solution
    differences = wd_grid - analytical_depth
    
    # Calculate L2 norm (RMS) = sqrt(sum(differences^2) / N)
    l2_norm = np.sqrt(np.sum(differences**2) / n_faces)
    
    grid_levels.append(grid_num + 1)
    l2_norms.append(l2_norm)
    grid_sizes.append(n_faces)
    
    print(f"Grid {grid_num+1}: {n_faces} faces, L2 norm = {l2_norm:.6e} m")

# Plot convergence
ax_convergence.loglog(grid_sizes, l2_norms, 'bo-', linewidth=2, markersize=8, label='L2 Norm')

# Add grid level annotations
for i, (size, norm, level) in enumerate(zip(grid_sizes, l2_norms, grid_levels)):
    ax_convergence.annotate(f'Grid {level}', 
                           xy=(size, norm), 
                           xytext=(10, 10), 
                           textcoords='offset points',
                           fontsize=10,
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))

# Add theoretical convergence lines for reference
if len(grid_sizes) >= 2:
    # Calculate slopes for 1st and 2nd order convergence
    x_ref = np.array([min(grid_sizes), max(grid_sizes)])
    
    # Reference lines (normalized to pass through first point)
    norm_ref = l2_norms[0] * (grid_sizes[0] / x_ref)**0.5  # 1st order (slope = -0.5)
    norm_ref2 = l2_norms[0] * (grid_sizes[0] / x_ref)**1.0  # 2nd order (slope = -1.0)
    
    ax_convergence.loglog(x_ref, norm_ref, 'r--', alpha=0.7, label='1st Order (slope = -0.5)')
    ax_convergence.loglog(x_ref, norm_ref2, 'g--', alpha=0.7, label='2nd Order (slope = -1.0)')

ax_convergence.set_xlabel('Number of Grid Faces')
ax_convergence.set_ylabel('L2 Norm (RMS Error) [m]')
ax_convergence.set_title('Grid Convergence Analysis: Water Depth Error vs Grid Resolution', fontweight='bold')
ax_convergence.grid(True, alpha=0.3, which='both')
ax_convergence.legend()

# Add text box with convergence rate calculation
if len(l2_norms) >= 2:
    # Calculate observed convergence rate between last two grids
    convergence_rate = np.log(l2_norms[-1] / l2_norms[-2]) / np.log(grid_sizes[-1] / grid_sizes[-2])
    
    textstr = f'Observed convergence rate\n(last two grids): {convergence_rate:.2f}'
    props = dict(boxstyle='round', facecolor='lightblue', alpha=0.7)
    ax_convergence.text(0.05, 0.95, textstr, transform=ax_convergence.transAxes, 
                       fontsize=10, verticalalignment='top', bbox=props)

plt.tight_layout()

# Save the convergence figure
plt.savefig(os.path.join(FIGURES, 'convergence_analysis.png'), dpi=300, bbox_inches='tight')
plt.show()

print(f"Convergence plot saved to: {os.path.join(FIGURES, 'convergence_analysis.png')}")
print("\n=== CONVERGENCE ANALYSIS ===")
for i, (level, size, norm) in enumerate(zip(grid_levels, grid_sizes, l2_norms)):
    print(f"Grid {level}: {size:4d} faces, L2 norm = {norm:.6e} m")
    if i > 0:
        rate = np.log(l2_norms[i] / l2_norms[i-1]) / np.log(grid_sizes[i] / grid_sizes[i-1])
        print(f"         Convergence rate from Grid {level-1}: {rate:.3f}")
print("==============================")
