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

# General parameters
g = 9.81        # gravitational acceleration (m/s^2)  
m_per_km = 1e3 # convert meters to kilometers

# For plotting the orders
xvgl1 = np.array([1e-3, 1e-2, 1e-1])
yvgl1 = xvgl1 ** (-1) / 1e5 # Minus 1 slope
yvgl2 = xvgl1 ** (-2) / 1e5 # Minus 2 slope

# Functions
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

def belanger_target_at_outflow(Q, C, L, B, ib, hout, d_outflow, g, dxfm=None, dc_grids=None, dx=0.5):
    """
    Compute the semi-analytical Belanger solution for water depth and level.
    
    Parameters:
    - Q: discharge (m³/s)
    - C: Chezy coefficient
    - L: channel length (m)
    - B: channel width (m)
    - ib: bed slope
    - hout: target water level at outflow (m)
    - d_outflow: target water depth at outflow (m)
    - g: gravitational acceleration (m/s²)
    - dxfm: cell size for single grid (optional)
    - dc_grids: list of cell sizes for multiple grids (optional)
    - dx: spatial step for Belanger integration (default 0.5 m)
    
    Returns:
    - If dc_grids is None: x, d, h arrays for the single grid
    - If dc_grids is provided: lists of x, d, h arrays for each grid
    """
    def compute_belanger(dxfm_local):
        dg = (Q**2 / B**2 / g)**(1/3)
        de = (Q**2 / B**2 / C**2 / ib)**(1/3)
        b_outflow = hout - d_outflow
        x = np.arange(0, L + dxfm + dx, dx)
        d = np.zeros_like(x)
        d[-1] = d_outflow
        for i in range(len(x) - 2, -1, -1):
            rhs = ib * (d[i + 1]**3 - de**3) / (d[i + 1]**3 - dg**3)
            d[i] = d[i + 1] - dx * rhs
        h = np.zeros_like(x)
        for i in range(len(x)):
            bed = b_outflow + ib * (x[-1] - x[i] - dxfm_local)
            h[i] = d[i] + bed
        return x, d, h

    if dc_grids is not None:
        xb_grids, deq_grids, heq_grids = [], [], []
        for dxfm_local in dc_grids:
            x, d, h = compute_belanger(dxfm_local)
            xb_grids.append(x)
            deq_grids.append(d)
            heq_grids.append(h)
        return xb_grids, deq_grids, heq_grids
    else:
        if dxfm is None:
            raise ValueError("Either dxfm or dc_grids must be provided.")
        return compute_belanger(dxfm)

def find_computation_points_in_analytical_solution(x_input, xb, deq, grids=None):
    """
    Finds the analytical solution at computation points.
    
    Parameters:
    - x_input: array of x-coordinates or full xcoord array if using grids
    - xb: array of x-coordinates where deq is defined
    - deq: array of analytical solution values corresponding to xb
    - grids: optional list of index arrays for subgrids
    
    Returns:
    - If grids is None: array of deq values corresponding to x_input
    - If grids is provided: list of arrays of deq values for each grid
    """
    def match_deq(x_nodes):
        deqs = np.zeros_like(x_nodes, dtype=float)
        for i, x in enumerate(x_nodes):
            j = np.where(xb == x)[0]
            if len(j) > 0:
                deqs[i] = deq[j[0]]
        return deqs

    if grids is None:
        return match_deq(x_input)
    else:
        deqs_grids = []
        for grid in grids:
            x_grid = x_input[grid]
            x_unique, _ = np.unique(x_grid, return_index=True)
            deqs_grids.append(match_deq(x_unique))
        return deqs_grids

def froude_number(Q, B, d):
    U = Q / (B * d)
    Fr = U / np.sqrt(g * d)
    return Fr

def compute_outflow_conditions(L, dxfm, ib):
    h_outflow = 0  # outflow water level w.r.t. reference
    b_outflow = -ib * (L + dxfm)  # bed level at outflow bc
    d_outflow = h_outflow - b_outflow  # water depth at outflow bc
    return h_outflow, b_outflow, d_outflow

def read_water_level_from_netcdf(dataset):
    xcoord_faces = dataset.variables['mesh2d_face_x'][:]  # Characteristic x-coordinate of mesh face
    ycoord_faces = dataset.variables['mesh2d_face_y'][:]  # Characteristic x-coordinate of mesh face
    # hfm = dataset.variables['mesh2d_s1'][:]
    hfm = dataset.variables['mesh2d_s1'][1][:] # Water level [m] 1 for selecting the latest time step
    return xcoord_faces, ycoord_faces, hfm

def compute_water_depth_fm_solution(xcoord_faces, dxfm, ib, hfm):
    xcoord_nodes = xcoord_faces + dxfm / 2  # Locations for water levels actually at netnodes
    ds = -xcoord_nodes * ib  # Compute water depth for fm solution
    dfm = hfm - ds
    return xcoord_nodes, dfm

def find_grids(ylink):
    grid1 = np.where((ylink >= 0) & (ylink <= 500))[0]
    grid2 = np.where((ylink >= 1000) & (ylink <= 1500))[0]
    grid3 = np.where((ylink >= 2000) & (ylink <= 2500))[0]
    grid4 = np.where((ylink >= 3000) & (ylink <= 3500))[0]
    grid5 = np.where((ylink >= 4000) & (ylink <= 4500))[0]
    return [grid1, grid2, grid3, grid4, grid5]

def calculate_rms_difference(deqs, dfm):
    ddiff = np.sqrt(np.mean((deqs - dfm) ** 2))
    return ddiff

def determine_l2_norm(paramfm, paramana, grids=None):
    """
    Computes the L2-norm of the difference between paramfm and paramana.
    
    Parameters:
    - paramfm: array or list of arrays (numerical solution)
    - paramana: array or list of arrays (analytical solution)
    - grids: optional list of index arrays if paramfm and paramana are full arrays
    
    Returns:
    - diffd: array of L2-norms for each grid
    """
    if grids is None:
        # Assume paramfm and paramana are lists of arrays
        diffd = np.zeros(len(paramfm))
        for i in range(len(paramfm)):
            diff = paramfm[i] - paramana[i]
            diffd[i] = np.sqrt(np.mean(diff ** 2))
    else:
        # Use grids to index into full arrays
        diffd = np.zeros(len(grids))
        for i, grid in enumerate(grids):
            diff = paramfm[grid] - paramana[grid]
            diffd[i] = np.sqrt(np.mean(diff ** 2))
    return diffd

def get_param_per_grid(param, grids):
    """
    Returns a list of arrays, each containing the values of 'param' for the indices in each grid.
    param: 1D array (e.g., dfm_faces, hfm, xcoord_faces, etc.)
    grids: list of index arrays
    """
    return [param[grid] for grid in grids]

def compute_water_depth_fm_solution_variable_dx(xcoord_faces, dc, ib, hfm):
    """
    Compute water depth at nodes for variable cell sizes.
    xcoord_faces: array of face x-coordinates (cell centers)
    dc: array of cell lengths (same length as xcoord_faces)
    ib: bed slope
    hfm: water level at cell centers
    Returns:
        xcoord_nodes: node x-coordinates (cell centers + 0.5*dc)
        dfm: water depth at nodes
    """
    xcoord_nodes = xcoord_faces + 0.5 * dc  # Node positions for variable cell sizes
    ds = -xcoord_nodes * ib
    dfm = hfm - ds
    return dfm

def average_and_center_interp(dfm_grids, ycoord_grids):
    """
    Computes the average of dfm_grids in the y-direction for each grid,
    and interpolates the result at the center of the domain in the y-direction.
    Returns:
        avg_per_grid: list of averages for each grid
        interp_center: interpolated value at y-center for each grid
    """
    avg_per_grid = []
    interp_center = []
    for dfm, y in zip(dfm_grids, ycoord_grids):
        avg_per_grid.append(np.mean(dfm))
        # Interpolate at y-center
        y_center = 0.5 * (np.min(y) + np.max(y))
        interp_val = np.interp(y_center, y, dfm)
        interp_center.append(interp_val)
    return avg_per_grid, interp_center

def average_dfm_y_direction_per_grid(xcoord_grids, ycoord_grids, dfm_grids,dc_grids=None):
    """
    For each grid, computes the average of dfm in the y-direction for each unique x-coordinate,
    and interpolates the value at the center of the y-domain for each x.
    Returns:
        x_unique_grids: list of arrays, each array is the unique x-coordinates for that grid
        avg_dfm_grids: list of arrays, each array is the y-averaged dfm for each unique x in that grid
        center_interp_grids: list of arrays, each array is the interpolated dfm at y-center for each unique x in that grid
    """
    avg_dfm_grids = []
    x_unique_grids = []
    center_interp_grids = []
    for xg, yg, dfmg in zip(xcoord_grids, ycoord_grids, dfm_grids):
        x_unique = np.unique(xg)
        avg_dfm = np.zeros_like(x_unique, dtype=float)
        center_interp = np.zeros_like(x_unique, dtype=float)
        y_center = 0.5 * (np.min(yg) + np.max(yg))
        for i, xval in enumerate(x_unique):
            mask = (xg == xval)
            avg_dfm[i] = np.mean(dfmg[mask])
            # Interpolate dfm at y_center for this x
            y_vals = yg[mask]
            dfm_vals = dfmg[mask]
            if len(y_vals) > 1:
                center_interp[i] = np.interp(y_center, y_vals, dfm_vals)
            elif len(y_vals) == 1:
                center_interp[i] = dfm_vals[0]
            else:
                center_interp[i] = np.nan
        avg_dfm_grids.append(avg_dfm)
        x_unique_grids.append(x_unique)
        center_interp_grids.append(center_interp)
    return x_unique_grids, avg_dfm_grids, center_interp_grids

def compute_observed_order_of_accuracy_finest(dc, diffd):
    """
    Computes the observed order of accuracy based on the given cell lengths and L2-norm differences.
    Returns:
        orders: array of observed orders of accuracy based on the finest grid.
        This is useful for comparing the order of accuracy between two different methods or parameters.
        dc: array of cell lengths
        diffd: array of L2-norm differences for each grid
        The order is computed relative to the finest grid (last element in dc). 
    """
    orders = np.zeros(len(dc) - 1)
    for i in range(len(dc) - 1):
        orders[i] = np.log(diffd[i] / diffd[-1]) / np.log(dc[i] / dc[-1])
    return orders

def plot_waterdepth_vs_belanger(
    xb, deq, xcoord, dfm, ddiff=None, loc_outflow_bc=None, d_outflow=None, L=None,
    method='', case_name='', is_grids=True, grid_labels=None, fig_folder='.', dxfm=0, dxbelanger=0.5, m_per_km = 1000.0
):
    """
    Plots water depth vs Belanger solution for one or multiple grids.
    If is_grids=True, xcoord and dfm should be lists of arrays (one per grid).
    If is_grids=False, xcoord and dfm should be single arrays.
    ddiff: list of rms differences for each grid or a single value.
    """
    plt.figure(1)
    plt.plot(xb / m_per_km, deq, 'k-', linewidth=lw, label='semi-analytical solution')

    if is_grids:
        colors = plt.cm.viridis(np.linspace(0, 1, len(xcoord)))
        for i in range(len(xcoord)):
            label = f'Grid {i+1} - {method}' if not grid_labels else grid_labels[i]
            plt.plot(xcoord[i] / m_per_km, dfm[i], '-', color=colors[i], linewidth=lw, markersize=ms2, label=label)
            if ddiff is not None:
                y_position = np.nanmean(dfm[i]) + i * 0.05 - 0.5
                plt.text(0.0, y_position, f'rms={ddiff[i]:.2e}', color=colors[i], fontsize=fs-4)
            y_min, y_max = np.min(dfm[i]), np.max(dfm[i])
    else:
        plt.plot(xcoord / m_per_km, dfm, 'r-', linewidth=lw, markersize=ms2, label=f'dflowfm, dx = {dxfm} m')
        if ddiff is not None:
            plt.text(16.5, 11.7, f'rms difference = {ddiff:.2e} m', color=[cl, cl, cl], fontsize=fs)
        y_min, y_max = np.min(dfm), np.max(dfm)
    if loc_outflow_bc is not None and d_outflow is not None:
        plt.plot(loc_outflow_bc / m_per_km, d_outflow, 'k.', linewidth=lw, markersize=ms3, label='boundary condition')

    plt.grid(True, linestyle=':', axis='both')
    plt.title('Water depth along the channel' + ('' if not is_grids else ' (per grid)'), fontsize=fs)
    if L is not None:
        plt.xlim([0, (L+dxfm*2)/m_per_km])
    plt.xlabel('distance from inflow boundary [km]', fontsize=fs)
    plt.gca().yaxis.set_major_locator(MultipleLocator(0.5))
    plt.gca().tick_params(axis='both', which='major', labelsize=(fs-3))

    plt.ylim(np.floor(y_min * 2) / 2, np.ceil(y_max * 2) / 2)
    plt.ylabel('water depth [m]', fontsize=fs)
    plt.legend(loc='best')
    plt.tight_layout()

    # Save the figure
    figure_name = f'{case_name}_water_depth_dxfm_{dxfm}_m_{method}_dxbelanger_{dxbelanger}_m'
    plt.savefig(f"{fig_folder}/{figure_name}.png", dpi=300)
    plt.close()
    
def plot_param_convergence(dc, diffd, param_name,method='',case_name=''):
    plt.figure(1)
    
    # Define the number of grids
    num_grids = len(dc)
    
    # Create a colormap
    cmap = plt.get_cmap('viridis', num_grids)
    
    # Plot each point with a different color
    for i in range(num_grids):
        plt.loglog(1.0 / dc[i], diffd[i], 'o', color=cmap(i), markersize=ms2, linewidth=lw, label=f'Grid {i+1} - {method}')
    
    # Plot the slopes
    plt.loglog(xvgl1, yvgl1, '--', color=[cl, cl, cl], linewidth=lw) # first order slope
    plt.loglog(xvgl1, yvgl2, '--', color=[cl, cl, cl], linewidth=lw) # second order slope


    # Add text annotations along the angle of the slopes
    plt.text(xvgl1[1], yvgl1[1], 'First Order Slope', fontsize=12, color='black', ha='right', rotation=-15)
    plt.text(xvgl1[1], yvgl2[1], 'Second Order Slope', fontsize=12, color='black', ha='right', rotation=-30.7)
    
    # Add grid, limits, labels, and title
    plt.grid(True, which = 'both', linestyle=':', axis='both')
    plt.xlim([1e-3, 1e-1])
    plt.ylim([1e-4, 1e1])
    plt.legend(loc='best')
    plt.xlabel(r'inverse cell length [m$^{-1}$]', fontsize=fs)
    plt.ylabel(r'L$_2$-norm [m]', fontsize=fs)
    plt.title(f'Convergence behavior ({param_name})', fontsize=fs)
    
    # Adjust layout and show plot
    plt.tight_layout()
    # plt.show()
    
    # Save the figure
    figure_name = f'{case_name}_{param_name}_fm_{method}_m_dxbelanger_{dxbelanger}_m'
    plt.savefig(fig_folder + '/' + figure_name + '.png', dpi=300)
    plt.close()
    
def plot_observed_order_of_accuracy(dc, orders1, orders2, param1, param2, method='', case_name=''):
    """
    Plots the observed order of accuracy.
    dc: array of cell lengths
    orders: array of observed orders of accuracy
    """
    plt.figure(figsize=(8, 6))
    plt.plot(dc[:-1], orders1, 'o', color='blue', markersize=ms2, linewidth=lw, label=f'{param1} - {method}')
    plt.plot(dc[:-1], orders2, 'd', color='orange', markersize=ms2, linewidth=lw, label=f'{param2} - {method}')
    
    # Add a horizontal line at y=1 for reference
    plt.axhline(y=1, color='red', linestyle=None, linewidth=lw, label='Expected order = 1')
    
    plt.xscale('log')
    plt.xlabel('Cell length (m)', fontsize=fs)
    plt.ylabel('Observed order of accuracy (-)', fontsize=fs)
    plt.title(f'{case_name}', fontsize=fs)
    plt.grid(True, linestyle=':', which='both')
    plt.legend(loc='best')
    
    # Adjust layout and show plot
    plt.tight_layout()
    # plt.show()
    
    # Save the figure
    figure_name = f'{case_name}_observed_order_of_accuracy_{method}'
    plt.savefig(fig_folder + '/' + figure_name + '.png', dpi=300)
    plt.close() 
    

if __name__=="__main__":
    # Define folders
    case_folder = 'C:/Users/star_kj/OneDrive - Stichting Deltares/Documents//04_Software/D-HYDRO/Validation_cases/'
    output_folder = 'dflowfmoutput'
    # File that constains the results to be analysed
    map_file = 'simplechannel_map.nc'  
    grid_file = 'simplechannel_net.nc'  
    
    """Section 3.1 """
    case_name = 'c010_belanger'
    fig_folder = case_folder + case_name + '/Figures'
    dataset = nc.Dataset(case_folder + case_name + '/' + output_folder + '/'+ map_file)
    
    # Input parameters for the channel - C010
    C = 65          # Chézy friction factor (m^0.5/s)
    Q = 600         # discharge (m^3/s)
    L = 100000      # length of the channel (m)
    B = 20          # channel width (m)
    ib = 1e-4       # bed slope
    dxfm = 500      # cell length (m)
    dxbelanger = 0.5 # discrete element length for Belanger equation (m)
    
    h_outflow, b_outflow, d_outflow = compute_outflow_conditions(L, dxfm, ib)
    xcoord_faces, ycoord_faces, hfm = read_water_level_from_netcdf(dataset)
    xcoord_nodes, dfm = compute_water_depth_fm_solution(xcoord_faces, dxfm, ib, hfm)
    
    # Compute analytical solution
    xb, deq, heq = belanger(Q, C, L, B, ib, d_outflow, g, dxfm, dxbelanger) #xb: x coordinates Belanger solution

    deqs = find_computation_points_in_analytical_solution(xcoord_nodes, xb, deq)
    ddiff = calculate_rms_difference(deqs, dfm)

    # Make plot waterdepth_vs_belanger
    loc_outflow_bc = L + dxfm # The location for the water level outflow bc is dx/2 outside the grid (mirrored location)
    
    plot_waterdepth_vs_belanger(
    xb, deq, xcoord_nodes, dfm, ddiff, loc_outflow_bc, d_outflow, L,
    method='', case_name=case_name, is_grids=False, grid_labels=None, fig_folder=fig_folder, dxfm=dxfm, dxbelanger=dxbelanger)

    
    
    """Section 3.2 """
    case_name = 'c012_Belanger_with_refinement_squares'
    fig_folder = case_folder + case_name + '/Figures'
    dataset = nc.Dataset(case_folder + case_name + '/' + output_folder + '/'+ map_file)
    gridset = nc.Dataset(case_folder + case_name + '/' + grid_file)
    # Input parameters for the channel - C012
    C = 65          # Chézy friction factor (m^0.5/s)
    Q = 2500        # discharge (m^3/s)
    L = 10000       # length of the channel (m)
    B = 500         # channel width (m)
    deltab = 1      # bed level difference (m)
    ib = deltab / L # bed slope
    loc_outflow_bc = L # The location for the water level outflow bc is dx/2 outside the grid (mirrored location)
    dxbelanger = 0.005 # discrete element length for Belanger equation (m)
    dout = 2.5  # chosen output water level
    dc = np.array([500.00, 250.00, 125.00, 62.50, 31.25])  # cell length (m)
    yloc_cross_section = 0.5 * B  # y-coordinate for cross-section
    
    # Read coordinates and water levels from the netCDF dataset
    xcoord_faces, ycoord_faces, hfm = read_water_level_from_netcdf(dataset)
    
    # Find the indices for each grid based on y-coordinates
    grids = find_grids(ycoord_faces) 
    
    # read bed_nodes
    bed_nodes = dataset.variables['mesh2d_node_z'][:]  # (n_nodes,) #node_z is the bed level at nodes
    # Get node indices for each face (edge)
    edge_nodes = dataset.variables['mesh2d_edge_nodes'][:]  # (n_edges, 2)

    
    # water levels at faces, but per grid
    hfm_grids = []
    xcoord_faces_grids = []
    ycoord_faces_grids = []
    for i, grid in enumerate(grids):
        dc_grid = np.full_like(xcoord_faces[grid], dc[i])  # Use the correct dc for this grid
        ycoord_faces_grids.append(ycoord_faces[grid])
        xcoord_faces_grids.append(xcoord_faces[grid])
        hfm_grids.append(hfm[grid])
    
    # Option 1: Read water depth at faces for each grid
    dfm_grids = []
    for i, grid in enumerate(grids):
        dc_grid = np.full_like(xcoord_faces[grid], dc[i])  # Use the correct dc for this grid
        dfm = dataset.variables['mesh2d_waterdepth'][:][1][grid]
        dfm_grids.append(dfm)
    
    # Option 2: Compute water depth at nodes for each grid
    dfm_grids = []
    for i, grid in enumerate(grids):
        dc_grid = np.full_like(xcoord_faces[grid], dc[i])  # Use the correct dc for this grid
        dfm = compute_water_depth_fm_solution_variable_dx(xcoord_faces[grid], dc_grid, ib, hfm[grid])
        dfm_grids.append(dfm)
    
    hout = 1.5      # target water level at outflow [m+NAP]
    d_outflow = 2.5 # target water depth at outflow [m]
    dxfm = 0.  # cell length for single grid (m)
    xb, deq, heq = belanger_target_at_outflow(Q, C, L, B, ib, hout, d_outflow, g, dxfm, dc_grids=None, dx=dxbelanger)
   
    # dc is your array of cell sizes per grid
    xb_grids, deq_grids, heq_grids = belanger_target_at_outflow(Q, C, L, B, ib, hout, d_outflow, g, dxfm, dc, dxbelanger)
    
    deqs_grids = []
    heqs_grids = []
    for i in range(len(grids)):
        deqs = find_computation_points_in_analytical_solution(
            xcoord_faces_grids[i], xb_grids[i], deq_grids[i], [np.arange(len(xcoord_faces_grids[i]))]
        )[0]
        heqs = find_computation_points_in_analytical_solution(
            xcoord_faces_grids[i], xb_grids[i], heq_grids[i], [np.arange(len(xcoord_faces_grids[i]))]
        )[0]
        deqs_grids.append(deqs)
        heqs_grids.append(heqs)
    
    # water depth
    x_unique_grids, avg_dfm_grids, center_interp_grids = average_dfm_y_direction_per_grid(xcoord_faces_grids, ycoord_faces_grids, dfm_grids, dc_grid)

    # Make plot (water depth convergence)
    method = 'average in y-direction'
    # ddiff_grids = [calculate_rms_difference(deqs_grids[i], avg_dfm_grids[i]) for i, grid in enumerate(grids)]
    ddiff_grids = determine_l2_norm(avg_dfm_grids, deqs_grids)
    plot_waterdepth_vs_belanger(xb, deq, x_unique_grids, avg_dfm_grids, ddiff_grids, loc_outflow_bc, d_outflow, L,
    method=method, case_name=case_name, is_grids=True, grid_labels=None, fig_folder=fig_folder, dxfm=dxfm, dxbelanger=dxbelanger)

    
    # diffd = np.array([calculate_rms_difference(avg_dfm_grids[i], deqs_grids[i]) for i in range(len(grids))])
    diffd = determine_l2_norm(avg_dfm_grids, deqs_grids)
    print(f"L2-norm for water depth: {diffd}")
    plot_param_convergence(dc, diffd, 'water depth',method,case_name)

    obs_order_of_acc_d_ave = compute_observed_order_of_accuracy_finest(dc, diffd)
    print(f"Observed order of accuracy: {obs_order_of_acc_d_ave}")
    
    method = f'cross-section at y = {round(yloc_cross_section)} m'
    # ddiff_grids = [calculate_rms_difference(deqs_grids[i], center_interp_grids[i]) for i, grid in enumerate(grids)]
    ddiff_grids = determine_l2_norm(center_interp_grids, deqs_grids)
    plot_waterdepth_vs_belanger(xb, deq, x_unique_grids, center_interp_grids, ddiff_grids, loc_outflow_bc, d_outflow, L,
    method=method, case_name=case_name, is_grids=True, grid_labels=None, fig_folder=fig_folder, dxfm=dxfm, dxbelanger=dxbelanger)
    
    # diffd = np.array([calculate_rms_difference(center_interp_grids[i], deqs_grids[i]) for i in range(len(grids))])
    diffd = determine_l2_norm(center_interp_grids, deqs_grids)
    print(f"L2-norm for water depth: {diffd}")
    plot_param_convergence(dc, diffd, 'water depth',method,case_name)
    
    obs_order_of_acc_d_center = compute_observed_order_of_accuracy_finest(dc, diffd)
    print(f"Observed order of accuracy: {obs_order_of_acc_d_center}")
    
    # Water level 
    x_unique_grids, avg_hfm_grids, center_interp_h_grids = average_dfm_y_direction_per_grid(xcoord_faces_grids, ycoord_faces_grids, hfm_grids)

    method = 'average in y-direction'
    # diffh = np.array([calculate_rms_difference(avg_hfm_grids[i], heqs_grids[i]) for i in range(len(grids))])
    diffh = determine_l2_norm(avg_hfm_grids, heqs_grids)
    print(f"L2-norm for water level: {diffh}")
    plot_param_convergence(dc, diffh, 'water level',method,case_name)
    
    obs_order_of_acc_h_ave = compute_observed_order_of_accuracy_finest(dc, diffh)
    print(f"Observed order of accuracy: {obs_order_of_acc_h_ave}")

    method = f'cross-section at y = {round(yloc_cross_section)} m'
    # diffh = np.array([calculate_rms_difference(center_interp_h_grids[i], heqs_grids[i]) for i in range(len(grids))])
    diffh = determine_l2_norm(center_interp_h_grids, heqs_grids)
    print(f"L2-norm for water level: {diffh}")
    plot_param_convergence(dc, diffh, 'water level',method,case_name)
    
    obs_order_of_acc_h_center = compute_observed_order_of_accuracy_finest(dc, diffh)
    print(f"Observed order of accuracy: {obs_order_of_acc_h_center}")

    method = f'cross-section at y = {round(yloc_cross_section)} m'
    plot_observed_order_of_accuracy(dc, obs_order_of_acc_d_center, obs_order_of_acc_h_center, 'water depth', 'water level', method, case_name)
    method = 'average in y-direction'
    plot_observed_order_of_accuracy(dc, obs_order_of_acc_d_ave, obs_order_of_acc_h_ave, 'water depth', 'water level', method, case_name)

