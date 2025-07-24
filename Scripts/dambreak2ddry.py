import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from PIL import Image
import netCDF4 as nc
import os


lw = 3
fs = 13
    
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
        xmin, xmax = 0, 26
        ymin, ymax = 0.3, 0.6
    elif meet == 1:
        xmin, xmax = 0, 26
        ymin, ymax = 0.0, 0.08
    elif meet == 2 or meet == 3 or meet == 4 or meet == 5 or meet == 6:
        xmin, xmax = 0, 26
        ymin, ymax = 0.0, 0.06
    else:
        xmin, xmax = 0, 26
        ymin, ymax = 0, 0.16


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
    ax.legend(loc='upper right', fontsize=fs-4, facecolor='white',framealpha=1)
    
    # ax.set_box_aspect(1)
    # plt.show()
    output_path = f'{figure_folder}/dambreak2ddry{meet:02d}.png'
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close(fig)
    
def plot_dambreak_contour(dataset_map, case_name, time_idx = 1, figure_folder='Figures'):   
    # Extract coordinates and water depth
    x = dataset_map.variables['mesh2d_face_x'][:]
    y = dataset_map.variables['mesh2d_face_y'][:]
    waterdepth = dataset_map.variables['mesh2d_waterdepth'][time_idx]  # Time index 1
    time_contour = np.round(dataset_map.variables['time'][time_idx])  # Time index 1
    # Plot using tricontourf for unstructured grid
    plt.figure(figsize=(10, 6))
    levels = np.linspace(0.0, 0.4, 10)
    contour = plt.tricontourf(x, y, waterdepth, levels=levels, cmap='viridis')
    plt.colorbar(contour, label='Water Depth (m)')
    # plt.clim([0, 0.3])
    plt.xlabel('X coordinate',fontsize=fs)
    plt.ylabel('Y coordinate',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.title(f'Water depth at time =  {time_contour} s',fontsize=fs)
    plt.grid(True, linestyle=':', axis='both')
    plt.tight_layout()

    
    output_path = f'{figure_folder}/dambreak2ddry_contour_time_{time_contour}_{case_name}.png'
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.show()



if __name__=="__main__":
    # Define folders
    case_folder = 'C:/Users/star_kj/OneDrive - Stichting Deltares/Documents/04_Software/D-HYDRO/Validation_cases/'
    output_folder = 'dflowfmoutput'
    his_file = 'dambreak2ddry_his.nc'  
    map_file = 'dambreak2ddry_map.nc'  
    case_name_squares = 'c031_dambreak2ddry_squares'
    case_name_triangles = 'c033_dambreak2ddry_triangles' #'c033_dambreak2ddry_triangles'
    fig_folder = case_folder + case_name_squares + '/Figures'

    dataset_squares = nc.Dataset(case_folder + case_name_squares + '/' + output_folder + '/'+ his_file)
    dataset_triangles = nc.Dataset(case_folder + case_name_triangles + '/' + output_folder + '/'+ his_file)
    
    dataset_squares_map = nc.Dataset(case_folder + case_name_squares + '/' + output_folder + '/'+ map_file)
    dataset_triangles_map = nc.Dataset(case_folder + case_name_triangles + '/' + output_folder + '/'+ map_file)
    
    if not os.path.exists(fig_folder):
        os.makedirs(fig_folder)     
    
    # Contour plots the dambreak simulation results for a specific time and domain of both cases
    plot_dambreak_contour(dataset_squares_map, case_name_squares, time_idx=1, figure_folder=fig_folder)
    plot_dambreak_contour(dataset_triangles_map, case_name_triangles, time_idx=1, figure_folder=fig_folder)

        
    # Plot the dambreak simulation results for a specific time and domain
    meet = 6 # location of the observation point
    im_path = f'{case_folder + case_name_squares}/dry/dbd_gauge{meet:02d}_meas.jpg'
    plot_dambreak(dataset_squares, dataset_triangles,im_path=im_path,  meet=meet, figure_folder=fig_folder)