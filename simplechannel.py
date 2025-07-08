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

    return x, d

# Froude number calculation
def froude_number(Q, B, d):
    U = Q / (B * d)
    Fr = U / np.sqrt(g * d)
    return Fr


# Input parameters for the channel
C = 65          # Chézy friction factor (m^0.5/s)
Q = 600         # discharge (m^3/s)
L = 100000      # length of the channel (m)
B = 20          # channel width (m)
ib = 1e-4       # bed slope
dxfm = 500      # cell length (m)
dxbelanger = 1000 # discrete element length for Belanger equation (m)
g = 9.81        # gravitational acceleration (m/s^2)  

# Boundary condition at outflow
h_outflow = 0 # outflow water level w.r.t.reference
b_outflow = - ib * (L + dxfm) # bed level at outflow bc
d_outflow = h_outflow - b_outflow # water depth at outflow bc

# Analytical solution
xb, deq = belanger(Q, C, L, B, ib, d_outflow, g, dxfm, dxbelanger) #xb: x coordinates Belanger solution


# Define folders
case_folder = 'C:/Users/star_kj/OneDrive - Stichting Deltares/Documents//04_Software/D-HYDRO/Validation_cases/c010_belanger'
output_folder = 'dflowfmoutput'
fig_folder = case_folder + '/Figures'

# File that constains the results to be analysed
map_file = 'simplechannel_map.nc'

# Read water level from NetCDF file
dataset = nc.Dataset(case_folder + '/' + output_folder + '/'+ map_file)
xcoord_faces = dataset.variables['mesh2d_face_x'][:] # Characteristic x-coordinate of mesh face
hfm = dataset.variables['mesh2d_s1'][:][-1, :]  #  standard_name: sea_surface_height, long_name: Water level [m]

# Locations for water levels actually at netnodes
xcoord_nodes = xcoord_faces + dxfm / 2

# Compute water depth for fm solution
ds = -xcoord_nodes * ib
dfm = hfm - ds

# Find computation points in analytical solution
deqs = np.zeros_like(xcoord_nodes)
for i in range(len(xcoord_nodes)):
    j = np.where(xb == xcoord_nodes[i])[0]
    if len(j) > 0:
        deqs[i] = deq[j[0]]

# Calculate RMS difference between semi-analytical and numerical solution
ddiff = np.sqrt(np.mean((deqs - dfm) ** 2))

# Calculate Froude number
Fr = froude_number(Q, B, deq)
print(f"Froude = {Fr}")

# Make plot
plt.figure(1)
plt.plot(xb / m_per_km, deq, 'k-', linewidth=lw, markersize=ms2, label=f'semi-analytical solution, dx = {dxbelanger} m')
plt.plot(xcoord_nodes / m_per_km, dfm, 'r-', linewidth=lw, markersize=ms2, label=f'dflowfm, dx = {dxfm} m')
loc_outflow_bc = L + dxfm # The location for the water level outflow bc is dx/2 outside the grid (mirrored location)
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
