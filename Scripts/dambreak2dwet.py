import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from PIL import Image
import netCDF4 as nc
import os

def plot_dambreak(dataset_squares, dataset_triangles, im_path, meet=3, figure_folder='Figures'):
    # Load data
    his = dataset_squares.variables['waterlevel'][:]
    time = dataset_squares.variables['time'][:]

    hist = dataset_triangles.variables['waterlevel'][:]
    timet = dataset_triangles.variables['time'][:]

    # Load image
    im = Image.open(im_path)

    # Dummy data
    dum = -1 * np.ones(len(time))

    # Axis limits
    if meet == 0:
        xmin, xmax = 0, 20
        ymin, ymax = 0.3, 0.6
    else:
        xmin, xmax = 0, 20
        ymin, ymax = 0, 0.16

    lw = 3
    fs = 13

    # Plot background image
    fig = plt.figure(figsize=(10, 6))
    plt.imshow(im)
    plt.axis('off')

    # Overlay plot
    ax = fig.add_axes([0.269, .195, .51, .631])
    ax.plot(time, his[:, meet], 'k', linewidth=lw, label='D-Flow FM squares')
    ax.plot(timet, hist[:, meet], 'm', linewidth=lw, label='D-Flow FM triangles')
    ax.plot(time, dum, 'r', linewidth=lw, label='Delft-FLS computation')
    ax.plot(time, dum, 'b', linewidth=lw, label='Stelling & Duinmeijer experiment (2001)')

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.set_xticks([xmin, xmax])
    ax.set_yticks([ymin, ymax])
    ax.set_facecolor('none')
    ax.legend(loc='lower right', fontsize=fs-1, facecolor='white',framealpha=1)
    
    # ax.set_box_aspect(1)
    # plt.show()
    output_path = f'{figure_folder}/dambreak2dwet{meet:02d}.png'
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close(fig)



if __name__=="__main__":
    # Define folders
    case_folder = 'C:/Users/star_kj/OneDrive - Stichting Deltares/Documents/04_Software/D-HYDRO/Validation_cases/'
    output_folder = 'dflowfmoutput'
    map_file = 'dambreak2dwet_his.nc'  
    case_name_squares = 'c030_dambreak2dwet_squares'
    case_name_triangles = 'c032_dambreak2dwet_triangles'
    fig_folder = case_folder + case_name_squares + '/Figures'

    dataset_squares = nc.Dataset(case_folder + case_name_squares + '/' + output_folder + '/'+ map_file)
    dataset_triangles = nc.Dataset(case_folder + case_name_triangles + '/' + output_folder + '/'+ map_file)
    
    # Create output folder if it doesn't exist
    if not os.path.exists(fig_folder):
        os.makedirs(fig_folder)     
        
    # Plot the dambreak simulation results for a specific time and domain
    meet = 0
    im_path = f'{case_folder + case_name_squares}/wet/dbw_gauge{meet:02d}_meas.jpg'
    plot_dambreak(dataset_squares, dataset_triangles,im_path=im_path,  meet=meet, figure_folder=fig_folder)