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
t_index = 0
t = times[t_index]

# Analytical solution for waterdepth(y)
T_sidereal = 23*3600 + 56*60 + 4.1  # seconds
Omega = 2 * np.pi / T_sidereal
phi = np.deg2rad(45)
f_coriolis = 2 * Omega * np.sin(phi)
y0 = 150000.0
g = 9.81

def compute_l2_norm(numerical, analytical):
    return np.sqrt(np.mean((numerical - analytical) ** 2))

def load_case(case_folder):
    dataset = nc.Dataset(os.path.join(case_folder, OUTPUT_FOLDER, MAP_FILE_NAME))
    data_vars = {var: dataset.variables[var][:] for var in dataset.variables}
    x = data_vars['mesh2d_face_x'][:]
    y = data_vars['mesh2d_face_y'][:]
    wd = data_vars['mesh2d_waterdepth'][:]
    u = data_vars['mesh2d_ucx'][:] if 'mesh2d_ucx' in data_vars else None
    v = data_vars['mesh2d_ucy'][:] if 'mesh2d_ucy' in data_vars else None
    bed_z = data_vars['mesh2d_face_z'][:] if 'mesh2d_face_z' in data_vars else np.zeros_like(x)
    times = data_vars['time'][:] if 'time' in data_vars else np.arange(wd.shape[0])
    return x, y, wd, u, v, bed_z, times

def get_cross_section(x, y, arr, x_center, tol, t_index):
    cross_idx = np.where(np.abs(x - x_center) < tol)[0]
    y_cross = y[cross_idx]
    arr_cross = arr[t_index, cross_idx]
    return y_cross, arr_cross, cross_idx

# Reference to all cases
base_folder = os.getcwd()
folders = [
    base_folder,
    os.path.join(base_folder, "..", "C030_flatbed_quadrilateral_refined"),
    os.path.join(base_folder, "..", "C030_flatbed_quadrilateral_doublerefined")
]
case_labels = ["Coarse", "Refined", "Double Refined"]

l2_waterlevel = []
l2_velocity = []

for folder in folders:
    x, y, wd, u, v, bed_z, times = load_case(folder)
    # Use first time step if only one is available
    t_index = 0

    # Water elevation for entire domain
    wd_domain = wd[t_index, :]
    bed_domain = bed_z
    elev_domain = wd_domain + bed_domain

    # Velocity magnitude for entire domain
    if u is not None and v is not None:
        u_domain = u[t_index, :]
        v_domain = v[t_index, :]
        vel_mag_domain = np.sqrt(u_domain**2 + v_domain**2)
        u_mean = np.mean(u_domain)
    elif u is not None:
        u_domain = u[t_index, :]
        vel_mag_domain = np.abs(u_domain)
        u_mean = np.mean(u_domain)
    else:
        vel_mag_domain = np.ones_like(wd_domain)
        u_mean = 1.0

    # Analytical solutions for entire domain
    T_sidereal = 23*3600 + 56*60 + 4.1
    Omega = 2 * np.pi / T_sidereal
    phi = np.deg2rad(45)
    f_coriolis = 2 * Omega * np.sin(phi)
    y0 = 150000.0
    g = 9.81

    wd_analytical_domain = -f_coriolis * (y - y0) / g * u_mean
    elev_analytical_domain = wd_analytical_domain + bed_domain
    vel_analytical_domain = np.full_like(vel_mag_domain, np.abs(u_mean))  # analytical velocity magnitude (constant)

    # Sort by y for norm calculation (optional, for consistency)
    sort_idx = np.argsort(y)
    elev_domain_sorted = elev_domain[sort_idx]
    elev_analytical_sorted = elev_analytical_domain[sort_idx]
    vel_mag_sorted = vel_mag_domain[sort_idx]
    vel_analytical_sorted = vel_analytical_domain[sort_idx]

    # L2 norms over the entire domain
    l2_waterlevel.append(compute_l2_norm(elev_domain_sorted, elev_analytical_sorted))
    l2_velocity.append(compute_l2_norm(vel_mag_sorted, vel_analytical_sorted))

# Example link lengths for each refinement
# Compute link lengths from grid data (structured cartesian grid)
link_lengths = []
for folder in folders:
    x, y, wd, u, v, bed_z, times = load_case(folder)
    # Assuming grid is structured along y, get unique sorted y values
    x_unique = np.unique(np.sort(x))
    # Link length is the minimum spacing between adjacent x values
    dx = np.min(np.diff(x_unique))
    link_lengths.append(dx)
link_lengths = np.array(link_lengths)
inv_link_lengths = 1 / link_lengths

print("Inverse link lengths:", inv_link_lengths)
# Plot convergence for water level L2 norm
plt.figure(figsize=(7, 5))
plt.loglog(inv_link_lengths, l2_waterlevel, 'o-', color='blue', label=r'$L_2$-norm (water level)')

# Theoretical first order line (choose a reference point)
ref_x = inv_link_lengths[0]
ref_y = l2_waterlevel[0]
first_order_y = ref_y * (inv_link_lengths**(-1)  / ref_x**(-1)) # slope = 1

plt.loglog(inv_link_lengths, first_order_y, 'k--', label='First order')

plt.xlabel('inverse flow link length [m$^{-1}$]')
plt.ylabel(r'$L_2$-norm [m]')
plt.title('Convergence behavior (water level)')
plt.grid(True, which='both', linestyle=':')
plt.legend()
plt.tight_layout()
plt.savefig("convergence_waterlevel.png", dpi=300, bbox_inches='tight')
plt.show()



