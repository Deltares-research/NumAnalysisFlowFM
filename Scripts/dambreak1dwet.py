import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import netCDF4 as nc
from matplotlib.ticker import MaxNLocator, MultipleLocator
import os

# Figure settings
fs = 14 # Font size

# General parameters
g = 9.81    # Gravitational acceleration in m/s^2

def propagationspeed(hl, hr):
    g = 9.81
    c1 = np.sqrt(g * hl)
    c0 = np.sqrt(g * hr)
    u20 = c0
    z0 = c0
    c20 = c0
    nter = 50
    for _ in range(nter):
        aa = 2 * z0 - u20
        ab = -c20
        ac = -z0
        ad = 0.5 * c0**2 + u20 * z0 - z0**2 + 0.5 * c20**2

        ba = -c20**2 + c0**2
        bb = 2 * c20 * u20 - 2 * c20 * z0
        bc = c20**2
        bd = -c20**2 * u20 + c20**2 * z0 - c0**2 * z0

        ca = 0.0
        cb = 2.0
        cc = 1.0
        cd = 2 * c1 - u20 - 2 * c20
        dd = aa * (bb * cc - bc * cb) - ab * (ba * cc - bc * ca) + ac * (ba * cb - bb * ca)
        d1 = ad * (bb * cc - bc * cb) - ab * (bd * cc - bc * cd) + ac * (bd * cb - bb * cd)
        d2 = aa * (bd * cc - bc * cd) - ad * (ba * cc - bc * ca) + ac * (ba * cd - bd * ca)
        d3 = aa * (bb * cd - bd * cb) - ab * (ba * cd - bd * ca) + ad * (ba * cb - bb * ca)

        dz = d1 / dd
        dc2 = d2 / dd
        du2 = d3 / dd
        z0 += dz
        c20 += dc2
        u20 += du2

    return c20


def dambreak1dwet(dataset, t=3600, domain=3):


    # Read data from NetCDF
    ncfile = dataset
    xdread = ncfile.variables['mesh2d_face_x'][:]
    ydread = ncfile.variables['mesh2d_face_y'][:]
    hdread = ncfile.variables['mesh2d_s1'][:]
    udread = ncfile.variables['mesh2d_ucx'][:]
    bdread = ncfile.variables['mesh2d_flowelem_ba'][:]

    # Define domain cases
    domain_bounds = {
        1: (-450, -250),
        2: (-250, -50),
        3: (-50, 150),
        4: (150, 450),
        5: (450, 650),
        6: (650, 950),
        7: (950, 5000)
    }
    yo, yb = domain_bounds.get(domain, (-500, 5000))

    # Select domain
    j = np.where((ydread > yo) & (ydread < yb))[0]
    xd = xdread[j] - L / 2
    yd = ydread[j]
    hd = hdread[:, j]
    ud = udread[:, j]
    bd = bdread[j]

    # Analytical solution
    cm = propagationspeed(hl, hr)
    hm = cm ** 2 / g

    xA = x0 - t * np.sqrt(g * hl)
    xB = x0 + 2 * t * np.sqrt(g * hl) - 3 * t * cm
    xC = x0 + t * (2 * cm ** 2 * (np.sqrt(g * hl) - cm)) / (cm ** 2 - g * hr)

    x1 = np.arange(0, xA, 1)
    h1 = hl * np.ones_like(x1)
    u1 = np.zeros_like(x1)

    x2 = np.arange(xA, xB, 1)
    h2 = 4 / 9 / g * (np.sqrt(g * hl) - (x2 - x0) / 2 / t) ** 2
    u2 = 2 / 3 * ((x2 - x0) / t + np.sqrt(g * hl))

    x3 = np.arange(xB, xC, 1)
    h3 = cm ** 2 / g * np.ones_like(x3)
    u3 = 2 * (np.sqrt(g * hl) - cm) * np.ones_like(x3)

    x4 = np.arange(xC, L, 1)
    h4 = hr * np.ones_like(x4)
    u4 = np.zeros_like(x4)

    x = np.concatenate([x1, x2, x3, x4]) - x0
    h = np.concatenate([h1, h2, h3, h4])
    u = np.concatenate([u1, u2, u3, u4])

    # Extract simulation data at time t
    time_index = -1
    hsom = hd[time_index, :]
    usom = ud[time_index, :]
    bsom = bd[:]

    # Sort data by x
    dat = np.column_stack((xd, yd, hsom, usom, bsom))
    dat = dat[np.argsort(dat[:, 0])]
    xd = dat[:, 0]
    yd = dat[:, 1]
    hsom = dat[:, 2]
    usom = dat[:, 3]
    bsom = dat[:, 4]

    # Plot water level
    plt.figure(1)
    plt.plot(x / 1000, h, 'k', linewidth=2, label='Analytical solution')
    plt.plot(xd / 1000, hsom, 'r', linewidth=2, label='D-Flow FM solution')
    plt.xlim([-30, 30])
    plt.ylim([0, 2.2])
    plt.grid(True)
    plt.legend(loc='lower left')
    plt.xlabel('Distance w.r.t. initial dam location [km]', fontsize=fs)
    plt.ylabel('Water level w.r.t. bed level [m]', fontsize=fs)
    plt.title(f'Surface level after t = {t:.0f} seconds', fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.tight_layout()
    figure_name = f'{case_name}_surface_level_t{int(t)}s_domain{domain}'
    plt.savefig(f"{fig_folder}/{figure_name}.png", dpi=300)
    plt.close()

    # Plot velocity
    plt.figure(2)
    plt.plot(x / 1000, u, 'k', linewidth=2, label='Analytical solution')
    plt.plot(xd / 1000, usom, 'r', linewidth=2, label='D-Flow FM solution')
    plt.xlim([-30, 30])
    plt.ylim([0, 5])
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.xlabel('Distance w.r.t. initial dam location [km]', fontsize=fs)
    plt.ylabel('Velocity in channel direction [m/s]', fontsize=fs)
    plt.title(f'Velocity in channel direction after t = {t:.0f} seconds', fontsize=fs)
    plt.xticks(fontsize=fs-2)
    plt.yticks(fontsize=fs-2)
    plt.tight_layout()
    figure_name = f'{case_name}_velocity_t{int(t)}s_domain{domain}'
    plt.savefig(f"{fig_folder}/{figure_name}.png", dpi=300)
    plt.close()

if __name__=="__main__":
    # Define folders
    case_folder = 'C:/Users/star_kj/OneDrive - Stichting Deltares/Documents/04_Software/D-HYDRO/Validation_cases/'
    output_folder = 'dflowfmoutput'
    map_file = 'dambreak1dwet_map.nc'  
    case_name = 'c020_dambreak1dwet_standard'
    fig_folder = case_folder + case_name + '/Figures'
    dataset = nc.Dataset(case_folder + case_name + '/' + output_folder + '/'+ map_file)

    # input parameters / Constants and initial conditions
    L = 60000   # Length of the channel in meters
    x0 = L / 2  # Initial dam location in meters

    hl = 2      # Left water level in meters
    hr = 0.1    # Right water level in meters
    
    # Create output folder if it doesn't exist
    if not os.path.exists(fig_folder):
        os.makedirs(fig_folder)     
        
    # Plot the dambreak simulation results for a specific time and domain
    dambreak1dwet(dataset, L, x0, hl, hr, t=3600.0,domain=6)