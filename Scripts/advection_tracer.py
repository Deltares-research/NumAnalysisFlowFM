import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset 
import netCDF4 as nc
import os

def read_nc_variable(filepath, varname):
    with Dataset(filepath, 'r') as nc:
        return np.array(nc.variables[varname][:])

def read_tracer_data(case_folder):
    data = {}
    data['Uadv'] = 0.1
    data['Trelease'] = 36000
    data['Tduration'] = 600
    data['printen'] = True
    # File paths
    paths = {
        'coarse': case_folder + '/c010_tracer_squares_coarse/dflowfmoutput/tracer_map.nc',
        'fine': case_folder + '/c011_tracer_squares_fine/dflowfmoutput/tracer_map.nc',
        'triangles': case_folder + '/c012_tracer_triangles/dflowfmoutput/tracer_map.nc'
    }

    # Read data
    for key, path in paths.items():
        data[f'block_{key}'] = read_nc_variable(path, 'mesh2d_block')
        data[f'time_{key}'] = read_nc_variable(path, 'time')
        data[f'cellcenter_{key}'] = read_nc_variable(path, 'mesh2d_face_x')

    return data

def plot_tracer(data, index=3851, figure_folder='Figures'):
    # index = 3851  corresponds to t = 2500 seconds after the release of the tracer
    fs = 16
    lw1 = 3
    lw2 = 3
    index2 = index #6417
    t = data['time_coarse'][index2]
    tracer1 = data['block_coarse'][index2, :]
    tracer2 = data['block_fine'][index, :]
    tracer3 = data['block_triangles'][index, :]

    fig = plt.figure(figsize=(10, 6))
    plt.plot(data['cellcenter_coarse'], tracer1, 'b', linewidth=lw1, label='coarse squares - one rows- Q')
    plt.plot(data['cellcenter_fine'], tracer2, 'g', linewidth=lw1, label='fine squares')
    plt.plot(data['cellcenter_triangles'], tracer3, 'r', linewidth=lw1, label='fine triangles')

    # Analytical solution
    xf = (t - data['Trelease']) * data['Uadv']
    xb = xf - data['Tduration'] * data['Uadv']
    plt.plot([-100, xb], [0, 0], 'k', linewidth=lw2)
    plt.plot([xb, xb], [0, 1], 'k', linewidth=lw2)
    plt.plot([xb, xf], [1, 1], 'k', linewidth=lw2)
    plt.plot([xf, xf], [1, 0], 'k', linewidth=lw2)
    plt.plot([xf, 500], [0, 0], 'k', linewidth=lw2, label='Analytical solution')

    # Beautify
    plt.xlim([-1, 300])
    plt.ylim([-0.1, 1.2])
    plt.grid(True, linestyle=':', axis='both')
    plt.xlabel('Distance from inflow boundary [m]', fontsize=fs)
    plt.ylabel('Tracer concentration [-]', fontsize=fs)
    time_after_release = t - data['Trelease']
    plt.title(f'Time is t = {time_after_release:.0f} seconds after release of the tracer', fontsize=fs)
    plt.legend(loc='upper left',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.tight_layout()
    plt.pause(0.02)

    output_path = f'{figure_folder}/tracer_time_{time_after_release:.0f}_coarse squares - one rows- Q.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)


if __name__ == "__main__":
    
     # Define folders
    case_folder = 'C:/Users/star_kj/OneDrive - Stichting Deltares/Documents/04_Software/D-HYDRO/Validation_cases/'
    output_folder = 'dflowfmoutput'
    map_file = 'tracer_map.nc'  
    case_name= 'c010_tracer_squares_coarse'
    fig_folder = case_folder + case_name + '/Figures'

    # Create output folder if it doesn't exist
    if not os.path.exists(fig_folder):
        os.makedirs(fig_folder)     
        
    # Plot the tracer simulation results
    data = read_tracer_data(case_folder)
    plot_tracer(data,index = 3850, figure_folder=fig_folder) #index = 3850
