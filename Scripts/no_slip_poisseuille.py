import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import netCDF4 as nc
import os


# Plot settings
fs = 16
ms = 20
    
def poiseuille_noslip_analysis(ncfile):
    # User-defined settings
    nuh = 0.1
    ny = 20
    B = 1e3
    L = 1e4
    deltaH = -1e-4

    # Physical constants and derived values
    g = 9.81
    deltay = B / ny
    dzdx = -deltaH / L

    # Analytical solution
    y = np.arange(-B / 2, B / 2 + deltay / 2, deltay / 2)
    uana = g * dzdx / (8 * nuh) * (B**2 - 4 * y**2)

    # Read simulation data from NetCDF
    # ncfile = Dataset('dflowfmoutput/poiseuillenoslip_his.nc', 'r')
    usim = ncfile.variables['x_velocity'][:][-1, :] # last time step
    vsim = ncfile.variables['y_velocity'][:][-1, :] # last time step
    # ncfile.close()

    # Compute differences
    diflon = (usim[ny // 2] - uana[ny]) / uana[ny] * 100
    diflat = np.max(np.abs(vsim))

    print(f'Difference near centerline is {diflon:.4f}%')
    print(f'Measure for the lateral velocities is {diflat:.4e}')


    # Plot analytical and simulated velocity profiles
    plt.figure()
    plt.plot(uana, y, label='Analytical solution', color='black', linewidth=2)
    plt.plot(usim, y[1::2], '.', markersize=ms, label='DFLOW-FM simulation', color='blue')
    plt.grid(True, linestyle=':', axis='both')
    plt.xlim([0.0, 0.15])
    plt.ylim([np.min(y), np.max(y)])
    plt.xlabel('Streamwise velocity [m/s]', fontsize=fs)
    plt.ylabel('Distance to channel centerline [m]', fontsize=fs)
    plt.yticks(np.arange(-B / 2, B / 2 + 1, 100))
    plt.legend(loc='best', fontsize=fs-4)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.title('Poiseuille flow with no-slip boundary condition', fontsize=fs)
    # plt.show()
    plt.tight_layout()
    figure_name = f'{case_name}_velocity'
    plt.savefig(f"{fig_folder}/{figure_name}.png", dpi=300)
    plt.close()
    
    
if __name__=="__main__":
    # Define folders
    case_folder = 'C:/Users/star_kj/OneDrive - Stichting Deltares/Documents//04_Software/D-HYDRO/Validation_cases/'
    output_folder = 'dflowfmoutput'
    # File that constains the results to be analysed
    map_file = 'poiseuillenoslip_map.nc'  
    his_file = 'poiseuillenoslip_his.nc'  
    
    case_name = 'c012_poiseuille_noslip'
    fig_folder = case_folder + case_name + '/Figures'
    
    # Create output folder if it doesn't exist
    if not os.path.exists(fig_folder):
        os.makedirs(fig_folder)     
        
    dataset = nc.Dataset(case_folder + case_name + '/' + output_folder + '/'+ his_file)

    poiseuille_noslip_analysis(dataset)
    