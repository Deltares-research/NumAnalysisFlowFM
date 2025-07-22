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


def dambreak1ddry(dataset, L, x0, hl, t=3600, domain=1):


    # Read data from NetCDF
    ncfile = dataset
    tdread = ncfile.variables['time'][:]
    xdread = ncfile.variables['mesh2d_face_x'][:]
    ydread = ncfile.variables['mesh2d_face_y'][:]
    hdread = ncfile.variables['mesh2d_s1'][:]
    udread = ncfile.variables['mesh2d_ucx'][:]
    # bdread = ncfile.variables['mesh2d_flowelem_ba'][:]

    # Define domain cases
    domain_bounds = {
        1: (-450, -250),
        2: (-250,  -50),
        3: ( -50,  150),
        4: ( 150,  450),
        5: ( 450,  650),
        6: ( 650,  950),
        7: ( 950, 5000)
    }
    yo, yb = domain_bounds.get(domain, (-500, 5000))

    # Pick case
    j = np.where((ydread > yo) & (ydread < yb))[0]
    xd = xdread[j] - L / 2
    yd = ydread[j]
    hd = hdread[:, j]
    ud = udread[:, j]

    # Divide into three domains
    xA = x0 - t * np.sqrt(g * hl)
    xB = x0 + 2 * t * np.sqrt(g * hl)

    # First part
    x1 = np.arange(0, xA + 1)
    h1 = hl * np.ones_like(x1)
    u1 = np.zeros_like(x1)

    # Second part: parabola
    x2 = np.arange(xA, xB + 1)
    h2 = 4 / 9 / g * (np.sqrt(g * hl) - (x2 - x0) / 2 / t) ** 2
    u2 = 2 / 3 * ((x2 - x0) / t + np.sqrt(g * hl))

    # Third part
    x3 = np.arange(xB, L + 1)
    h3 = np.zeros_like(x3)
    u3 = np.zeros_like(x3)

    # Collect data
    x = np.concatenate([x1, x2, x3]) - x0
    h = np.concatenate([h1, h2, h3])
    u = np.concatenate([u1, u2, u3])

    # Read data at specific time
    # time_idx = -1
    time_idx = np.where(tdread == t)[0][0]
    hsom = hd[time_idx, :]
    usom = ud[time_idx, :]

    # Sort data
    dat = np.column_stack((xd, yd, hsom, usom))
    dat = dat[np.argsort(dat[:, 0])]
    xd = dat[:, 0]
    yd = dat[:, 1]
    hsom = dat[:, 2]
    usom = dat[:, 3]

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
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.tight_layout()
    figure_name = f'{case_name}_water_level_t{int(t)}s_domain{domain}'
    plt.savefig(f"{fig_folder}/{figure_name}.png", dpi=300)
    plt.close()

    # Plot longitudinal velocity
    plt.figure(2)
    plt.plot(x / 1000, u, 'k', linewidth=2, label='Analytical solution')
    plt.plot(xd / 1000, usom, 'r', linewidth=2, label='D-Flow FM solution')
    plt.xlim([-30, 30])
    plt.ylim([0, 10])
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.xlabel('Distance w.r.t. initial dam location [km]', fontsize=fs)
    plt.ylabel('Velocity in channel direction [m/s]', fontsize=fs)
    plt.title(f'Velocity in channel direction after t = {t:.0f} seconds', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.tight_layout()
    figure_name = f'{case_name}_velocity_t{int(t)}s_domain{domain}'
    plt.savefig(f"{fig_folder}/{figure_name}.png", dpi=300)
    plt.close()


if __name__=="__main__":
    # Define folders
    case_folder = 'C:/Users/star_kj/OneDrive - Stichting Deltares/Documents/04_Software/D-HYDRO/Validation_cases/'
    output_folder = 'dflowfmoutput'
    map_file = 'dambreak1ddry_map.nc'  
    case_name = 'c021_dambreak1ddry_standard'
    fig_folder = case_folder + case_name + '/Figures'
    dataset = nc.Dataset(case_folder + case_name + '/' + output_folder + '/'+ map_file)

    # input parameters / Constants and initial conditions
    L = 60000   # Length of the channel in meters
    x0 = L / 2  # Initial dam location in meters

    hl = 2      # Left water level in meters
    hr = 0.0    # Right water level in meters
    
    # Create output folder if it doesn't exist
    if not os.path.exists(fig_folder):
        os.makedirs(fig_folder)     
        
    # Plot the dambreak simulation results for a specific time and domain
    dambreak1ddry(dataset, L, x0, hl,  t=3200.0,domain=6)