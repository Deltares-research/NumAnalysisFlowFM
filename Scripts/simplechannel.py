"""This script is a python script of the simplechannel.m script of test case c010_belanger"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from matplotlib.ticker import MaxNLocator, MultipleLocator

# Figure settings
ms1 = 15    # markersize1
ms2 = 8     # markersize2
ms3 = 12    # markersize3
fs = 16     # fontsize
cl = 0.5    # color
lw = 1      # linewidth
m_per_km = 1e3 # convert meters to kilometers

# General parameters
g = 9.81        # gravitational acceleration (m/s^2)  

# For plotting the orders
# Minus 1 slope
xvgl1 = np.array([1e-4, 1e0])
yvgl1 = xvgl1 ** (-1) / 1e5

# Minus 2 slope
xvgl2 = np.array([1e-4, 1e0])
yvgl2 = xvgl2 ** (-2) / 1e5

def belanger(Q, C, L, B, ib, hb, g, dxfm, dx = 0.5):
    """Compute the semi-analytical Belanger equation for *d* as the water depth 
    according to dflowfm_validation_doc.pdf Section 3.1"""

    # Compute the equilibrium depth *de* and the limit depth (associated with Fr = 1) *dg*
    dg = (Q**2 / B**2 / g)**(1/3)
    de = (Q**2 / B**2 / C**2 / ib)**(1/3)

    x = np.arange(0, L + dxfm + dx, dx)
    
    # Compute water depth d
    d = np.zeros(len(x))
    d[-1] = hb
    
    for i in range(len(x) - 2, -1, -1):
        rhs = ib * (d[i + 1]**3 - de**3) / (d[i + 1]**3 - dg**3)
        d[i] = d[i + 1] - dx * rhs
    
    h = np.zeros(len(x))
    for i in range(len(x)):
        h[i] = d[i] - ib * x[i]
    return x, d, h

def belanger_target_at_outflow(Q, C, L, B, ib, hout, d_outflow, g, dxfm, dx=0.5):
    """
    Compute the semi-analytical Belanger equation for *d* as the water depth.
    hout: target water level at outflow boundary (m)
    d_outflow: target water depth at outflow boundary (m)
    """
    # Compute the equilibrium depth *de* and the limit depth (associated with Fr = 1) *dg*
    dg = (Q**2 / B**2 / g)**(1/3)
    de = (Q**2 / B**2 / C**2 / ib)**(1/3)

    # Compute bed level at outflow
    b_outflow = hout - d_outflow

    x = np.arange(0, L + dxfm + dx, dx)

    # Compute water depth d
    d = np.zeros(len(x))
    # Set water depth at outflow boundary (last point)
    d[-1] = d_outflow

    for i in range(len(x) - 2, -1, -1):
        rhs = ib * (d[i + 1]**3 - de**3) / (d[i + 1]**3 - dg**3)
        d[i] = d[i + 1] - dx * rhs

    # Compute water level h = d + bed level at x
    h = np.zeros(len(x))
    for i in range(len(x)):
        bed = b_outflow + ib * (x[-1] - x[i])  # bed rises upstream
        h[i] = d[i] + bed
    return x, d, h


def compute_bedlevel_at_flowlink(dataset, bedlevtype=3):
    """
    Compute the bed level at each velocity point (cell face) as used in Delft3D-FM.
    bedlevtype:
        3: average of the two adjacent node elevations
        4: minimum of the two adjacent node elevations
        5: maximum of the two adjacent node elevations
    Returns:
        zk: array of bed levels at each face (np.nan for invalid edges)
        node_pairs: list of (node1, node2) indices for each face
    """
    edge_nodes = dataset.variables['mesh2d_edge_nodes'][:]  # (n_edges, 2)
    node_z = dataset.variables['mesh2d_node_z'][:]          # (n_nodes,)

    zk = np.full(edge_nodes.shape[0], np.nan)
    node_pairs = []

    for i, (n1, n2) in enumerate(edge_nodes):
        node_pairs.append((n1, n2))
        # Check for out-of-bounds or dummy indices
        if (n1 < 0 or n2 < 0 or n1 >= len(node_z) or n2 >= len(node_z)):
            continue
        if bedlevtype == 3:
            zk[i] = 0.5 * (node_z[n1] + node_z[n2])
        elif bedlevtype == 4:
            zk[i] = min(node_z[n1], node_z[n2])
        elif bedlevtype == 5:
            zk[i] = max(node_z[n1], node_z[n2])
        else:
            raise ValueError("Invalid bedlevtype. Use 3 (mean), 4 (min), or 5 (max).")
    return zk, node_pairs
       

# Froude number calculation
def froude_number(Q, B, d):
    U = Q / (B * d)
    Fr = U / np.sqrt(g * d)
    return Fr

# Make plot waterdepth_vs_belanger
def plot_waterdepth_vs_belanger(xb, deq, xcoord_nodes, dfm, loc_outflow_bc, d_outflow):
    plt.figure(1)
    plt.plot(xb / m_per_km, deq, 'k-', linewidth=lw, markersize=ms2, label=f'semi-analytical solution, dx = {dxbelanger} m')
    plt.plot(xcoord_nodes / m_per_km, dfm, 'r-', linewidth=lw, markersize=ms2, label=f'dflowfm, dx = {dxfm} m')
    plt.plot(loc_outflow_bc / m_per_km, d_outflow, 'k.', linewidth=lw, markersize=ms3, label='boundary condition')
    plt.grid(True,linestyle=':',axis='both')

    # Set labels and legend
    plt.title('Water depth along the channel', fontsize=fs)
    plt.xlim([0, 102])
    plt.xlabel('distance from inflow boundary [km]', fontsize=fs)
    # Set the y-axis to have ticks rounded to the nearest 0.5
    plt.gca().yaxis.set_major_locator(MultipleLocator(0.5))
    # Set the y-axis limits to be rounded at 0.5
    y_min, y_max = np.min(dfm), np.max(dfm)
    plt.ylim(np.floor(y_min * 2) / 2, np.ceil(y_max * 2) / 2)
    plt.ylabel('water depth [m]', fontsize=fs)
    plt.text(16.5, 11.7, f'rms difference = {ddiff:.2e} m', color=[cl, cl, cl], fontsize=fs)
    plt.gca().tick_params(axis='both', which='major', labelsize=fs)
    plt.legend(loc='best')
    plt.tight_layout()

    # Save plot
    figure_name = f'c010_belanger_waterdepth_dxfm_{dxfm}_m_dxbelanger_{dxbelanger}_m'
    plt.savefig(fig_folder + '/' + figure_name + '.png', dpi=300)
    plt.close()

# Make plot (param convergence)
def plot_param_convergence(dc, diffd, param_name):
    plt.figure(1)
    plt.loglog(1.0 / dc, diffd, 'o', markersize=ms2, linewidth=lw)
    plt.loglog(xvgl1, yvgl1, '--', color=[cl, cl, cl], linewidth=lw) # first order slope
    plt.loglog(xvgl2, yvgl2, '--', color=[cl, cl, cl], linewidth=lw) # second order slope
    plt.grid(True)
    plt.xlim([1e-3, 1e-1])
    plt.ylim([1e-4, 1e3])
    plt.legend(['Cartesian grids'], loc='best')
    plt.xlabel(r'inverse flow link length [m$^{-1}$]', fontsize=fs)
    plt.ylabel(r'L$_2$-norm [m]', fontsize=fs)
    plt.title(f'Convergence behavior ({param_name})', fontsize=fs)
    # plt.gca().set_fontsize(fs)
    plt.tight_layout()
    plt.show()
    figure_name = f'c010_belanger_{param_name}_dxfm_{dxfm}_m_dxbelanger_{dxbelanger}_m'
    plt.savefig(fig_folder + '/' + figure_name + '.png', dpi=300)
    plt.close()

def compute_outflow_conditions(L, dxfm, ib):
    h_outflow = 0  # outflow water level w.r.t. reference
    b_outflow = -ib * (L + dxfm)  # bed level at outflow bc
    d_outflow = h_outflow - b_outflow  # water depth at outflow bc
    return h_outflow, b_outflow, d_outflow

def read_water_level_from_netcdf(dataset):
    xcoord_faces = dataset.variables['mesh2d_face_x'][:]  # Characteristic x-coordinate of mesh face
    ycoord_faces = dataset.variables['mesh2d_face_y'][:]  # Characteristic x-coordinate of mesh face
    hfm = dataset.variables['mesh2d_s1'][:][-1, :]  # Water level [m]
    return xcoord_faces, ycoord_faces, hfm


def compute_water_depth_fm_solution(xcoord_faces, dxfm, ib, hfm):
    xcoord_nodes = xcoord_faces + dxfm / 2  # Locations for water levels actually at netnodes
    ds = -xcoord_nodes * ib  # Compute water depth for fm solution
    dfm = hfm - ds
    return xcoord_nodes, dfm

def find_computation_points_in_analytical_solution(xcoord_nodes, xb, deq):
    deqs = np.zeros_like(xcoord_nodes)
    for i in range(len(xcoord_nodes)):
        j = np.where(xb == xcoord_nodes[i])[0]
        if len(j) > 0:
            deqs[i] = deq[j[0]]
    return deqs

def find_computation_points_in_analytical_solution_grids(xcoord_nodes, xb, deq, grids):
    """
    For each grid, find the analytical solution at the computation points.
    Returns a list of arrays, one for each grid.
    """
    deqs_grids = []
    for grid in grids:
        deqs = np.zeros_like(grid, dtype=float)
        for idx, node_idx in enumerate(grid):
            x_node = xcoord_nodes[node_idx]
            j = np.where(xb == x_node)[0]
            if len(j) > 0:
                deqs[idx] = deq[j[0]]
        deqs_grids.append(deqs)
    return deqs_grids

def calculate_rms_difference(deqs, dfm):
    ddiff = np.sqrt(np.mean((deqs - dfm) ** 2))
    return ddiff


def find_grids(ylink):
    grid1 = np.where((ylink >= 0) & (ylink <= 500))[0]
    grid2 = np.where((ylink >= 1000) & (ylink <= 1500))[0]
    grid3 = np.where((ylink >= 2000) & (ylink <= 2500))[0]
    grid4 = np.where((ylink >= 3000) & (ylink <= 3500))[0]
    grid5 = np.where((ylink >= 4000) & (ylink <= 4500))[0]
    return [grid1, grid2, grid3, grid4, grid5]

def determine_l2_norm(paramfm, paramana, grids):
    diffd = np.zeros(len(grids))
    for i, grid in enumerate(grids):
        dhsqr = (paramfm[grid] - paramana[grid]) ** 2
        dhsqrsum = np.sum(dhsqr)
        diffd[i] = np.sqrt(dhsqrsum) / np.sqrt(len(grid))
    return diffd

def compute_water_depth_at_faces(hfm, zk, edge_nodes):
    """
    Compute water depth at faces using upstream water level minus bed level at face.
    hfm: water level at cell centers (faces are between centers)
    zk: bed level at faces
    edge_nodes: (n_edges, 2) array with node indices for each face
    Returns:
        dfm_faces: water depth at each face
    """
    dfm_faces = np.full(len(zk), np.nan)
    for i, (n1, n2) in enumerate(edge_nodes):
        # Choose upstream node (for 1D, the lower index is usually upstream)
        upstream = min(n1, n2)
        if upstream < 0 or upstream >= len(hfm):
            continue
        dfm_faces[i] = hfm[upstream] - zk[i]
    return dfm_faces

if __name__=="__main__":
    # Define folders
    case_folder = 'C:/Users/star_kj/OneDrive - Stichting Deltares/Documents//04_Software/D-HYDRO/Validation_cases/'
    output_folder = 'dflowfmoutput'
    # File that constains the results to be analysed
    map_file = 'simplechannel_map.nc'  
    
    # """Section 3.1 """
    # case_name = 'c010_belanger'
    # fig_folder = case_folder + case_name + '/Figures'
    # dataset = nc.Dataset(case_folder + case_name + '/' + output_folder + '/'+ map_file)
        
    # # Input parameters for the channel - C010
    # C = 65          # Chézy friction factor (m^0.5/s)
    # Q = 600         # discharge (m^3/s)
    # L = 100000      # length of the channel (m)
    # B = 20          # channel width (m)
    # ib = 1e-4       # bed slope
    # dxfm = 500      # cell length (m)
    # dxbelanger = 0.5 # discrete element length for Belanger equation (m)
    
    
    # h_outflow, b_outflow, d_outflow = compute_outflow_conditions(L, dxfm, ib)
    # xcoord_faces, ycoord_faces, hfm = read_water_level_from_netcdf(dataset)
    # xcoord_nodes, dfm = compute_water_depth_fm_solution(xcoord_faces, dxfm, ib, hfm)
    
    # # Compute analytical solution
    # xb, deq, heq = belanger(Q, C, L, B, ib, d_outflow, g, dxfm, dxbelanger) #xb: x coordinates Belanger solution

    # deqs = find_computation_points_in_analytical_solution(xcoord_nodes, xb, deq)
    # ddiff = calculate_rms_difference(deqs, dfm)

    # # Make plot waterdepth_vs_belanger
    # loc_outflow_bc = L + dxfm # The location for the water level outflow bc is dx/2 outside the grid (mirrored location)
    # plot_waterdepth_vs_belanger(xb, deq, xcoord_nodes, dfm, loc_outflow_bc, d_outflow)
    
    
    
    """Section 3.2 """
    case_name = 'c012_Belanger_with_refinement_squares'
    fig_folder = case_folder + case_name + '/Figures'
    dataset = nc.Dataset(case_folder + case_name + '/' + output_folder + '/'+ map_file)

    # Input parameters for the channel - C012
    C = 65          # Chézy friction factor (m^0.5/s)
    Q = 2500        # discharge (m^3/s)
    L = 10000       # length of the channel (m)
    B = 500         # channel width (m)
    deltab = 1
    ib = deltab / L # bed slope
    dxfm = 500      # cell length (m)
    dxbelanger = 0.1 # discrete element length for Belanger equation (m)
    dout = 2.5  # chosen output water level
    
    xcoord_faces, ycoord_faces, hfm = read_water_level_from_netcdf(dataset)
    
    
    # Compute l2 norm water depth
    xcoord_faces, ycoord_faces, hfm = read_water_level_from_netcdf(dataset)
    xcoord_nodes, dfm = compute_water_depth_fm_solution(xcoord_faces, dxfm, ib, hfm)
    h_outflow, b_outflow, d_outflow = compute_outflow_conditions(L, dxfm, ib)
    # xb, deq, heq = belanger(Q, C, L, B, ib, d_outflow, g, dxfm, dxbelanger) #xb: x coordinates Belanger solution
    
    hout = 1.5      # target water level at outflow [m+NAP]
    d_outflow = 2.5 # target water depth at outflow [m]
    xb, deq, heq = belanger_target_at_outflow(Q, C, L, B, ib, hout, d_outflow, g, dxfm, dxbelanger)
    
    deqs = find_computation_points_in_analytical_solution(xcoord_nodes, xb, deq)
    heqs = find_computation_points_in_analytical_solution(xcoord_nodes, xb, heq)
    
    zk, node_pairs = compute_bedlevel_at_flowlink(dataset, bedlevtype=3)
    edge_nodes = np.array(node_pairs)
    dfm_faces = compute_water_depth_at_faces(hfm, zk, edge_nodes)
    grids = find_grids(ycoord_faces)
    

    diffd = determine_l2_norm(dfm_faces, deqs, grids)
    diffh = determine_l2_norm(hfm, heqs, grids)
    print(f"L2-norm for water depth: {diffd}")

    # Cartesian grids
    dc = np.array([500.00, 250.00, 125.00, 62.50, 31.25])

    # # Make plot (water depth convergence)
    plot_param_convergence(dc, diffd, 'water depth')

    # # Make plot (velocity convergence)
    plot_param_convergence(dc, diffh, 'water level')

