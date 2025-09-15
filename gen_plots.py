#%%
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import sys
import os
sys.path.insert(0,os.getcwd())
from config import *
import glob
import libs.taper_and_rounding_functions as tf
import libs.field_functions as ff
import libs.array_functions as af

x_pts = np.linspace(-n2f_width/2,n2f_width/2,n_angles)
z_pts = np.linspace(-n2f_width/2,n2f_width/2,n_angles)

alphas = [np.rad2deg(np.arctan(f)) for f in x_pts/n2f_height]
betas = [np.rad2deg(np.arctan(f)) for f in z_pts/n2f_height]
aa,bb = np.meshgrid(alphas,betas)

#import wavelengths
temp_arr = []   #temp array to store the individial wavelength arrays
for file in glob.glob(f"{data_input_directory}/wvls-*.npy"):
    temp_arr.append(np.load(file))
all_wvls = np.concatenate(temp_arr,0)   #large array containing wavelengths from all sims sequentially

#import inc angles
temp_arr = []
for file in glob.glob(f"{data_input_directory}/inc_angles-*.npy"):
    temp_arr.append(np.load(file))
all_inc_angles = np.concatenate(temp_arr,0)     #large array containing incident angles from all sims sequentially

#import reflection data
#allow_pickle=True is required to import arrays/dictionaries containing objects themselves
#.item() is required because importing dictionaries doesn't work without it ¯\_(ツ)_/¯
temp_arr = []
if polarisation_setting != "mixed":
    for file in glob.glob(f"{data_input_directory}/refl-*.npy"):
        loaded = np.load(file,allow_pickle=True).item()
        pre_norm = np.absolute(loaded[ff.pol_string(polarisation_setting)])**2
        temp_arr.append(pre_norm)
else:
    for file in glob.glob(f"{data_input_directory}/refl-*.npy"):
        loaded = np.load(file,allow_pickle=True).item()
        pre_norm = np.absolute(loaded["Ex"])**2 + np.absolute(loaded["Ez"])**2
        temp_arr.append(pre_norm)
all_refl_coeffs = np.concatenate(temp_arr,2)     #large array containing reflection coefficients from all sims sequentially
print(all_refl_coeffs.shape)

#a working but unused script to take slices of the data corresponding to different incident angles. Will need a LOT of data to look halfway decent.
#also needs the commented parameters to be uncommented in ./config.py
'''
from mayavi import mlab
#plot slices corresponding to difference incident angles
inc_slice_vals = np.arange(inc_slice_min,inc_slice_max+inc_slice_interval,inc_slice_interval)
ax_ranges = [min(alphas),max(alphas),min(betas),max(betas),0,1]
ax_scale = [1.0,1.0,0.4]
ax_extent = ax_ranges * np.repeat(ax_scale,2)
for slice in inc_slice_vals:    #each loop is a new surface
    indices = af.get_indices_within(all_inc_angles,slice,0.5)   #each element is an array containing one element
    print(f"number of data points in slice with incident angle {slice}:")
    print(len(indices))
    slice_refl_coeffs = np.zeros((len(alphas),len(betas),len(indices)))
    slice_wvls = np.zeros((len(indices)))
    interp_peak_wvls = np.zeros((len(alphas),len(betas)))   #for this slice only
    interp_max_Rs = np.zeros((len(alphas),len(betas)))      
    for i in range(len(indices)):   #take slice at the chosen incident angle
        slice_refl_coeffs[:,:,i] = all_refl_coeffs[:,:,indices[i][0]]
        slice_wvls[i] = all_wvls[indices[i][0]]
    for i in range(slice_refl_coeffs.shape[0]):
        for j in range(slice_refl_coeffs.shape[1]):
            max3indices = np.argpartition(slice_refl_coeffs[i,j,:],-3)[-3:]
            print(max3indices)
            R_maxes = (slice_refl_coeffs[i,j,max3indices[0]],slice_refl_coeffs[i,j,max3indices[1]],slice_refl_coeffs[i,j,max3indices[2]])
            print(R_maxes)
            lambda_maxes = (slice_wvls[max3indices[0]],slice_wvls[max3indices[1]],slice_wvls[max3indices[2]])
            print(lambda_maxes)
            interp = af.parabolic_max(lambda_maxes,R_maxes)
            interp_peak_wvls[i,j] = interp[0]
            interp_max_Rs[i,j] = interp[1]
    fig = mlab.figure()
    mlab.surf(alphas,betas,interp_max_Rs,warp_scale="auto")   
    mlab.savefig(filename=f"int-surface-{slice}.png") 
    mlab.show()
    fig = mlab.figure()
    mlab.surf(alphas,betas,interp_peak_wvls,warp_scale="auto")
    mlab.savefig(filename=f"wvl-surface-{slice}.png")
    mlab.show()
print("Done with slices!")
'''

###ensemble average across all the incident angles used (for each wvl)

avgs = []
avg_wvls = np.zeros(n_freq)
for i in range(n_freq):
    this_wvl_inds = np.where(all_wvls == all_wvls[i])[0]
    this_wvl_avg = np.zeros((len(alphas),len(betas),1))
    for k in this_wvl_inds:
        this_wvl_avg += np.atleast_3d(all_refl_coeffs[:,:,k])
    this_wvl_avg = np.divide(this_wvl_avg,len(this_wvl_inds))
    avgs.append(this_wvl_avg)
    avg_wvls[i] = all_wvls[i]
avg_refl_coeffs = np.concatenate(avgs,2)    #data across all bloch angles averaged into one array
print("successfully made ensemble average")
print(f"number of wavelengths: {avg_refl_coeffs.shape[2]}")

integrated_Rs = np.zeros((len(alphas),len(betas)))  #total integrated power (wrt wavelength) for each pair of reflected angles
for i in range(avg_refl_coeffs.shape[0]):
    for j in range(avg_refl_coeffs.shape[1]):
        integral = np.trapezoid(avg_refl_coeffs[i,j,:],avg_wvls)
        integrated_Rs[j,i] = integral

interp_peak_wvls = np.zeros((len(alphas),len(betas)))   #parabolically interpolated wavelength (corresponding to max power) for each pair of reflected angles.
interp_max_Rs = np.zeros((len(alphas),len(betas)))      #parabolically interpolated maximum power for each pair of reflected angles
for i in range(avg_refl_coeffs.shape[0]):
    for j in range(avg_refl_coeffs.shape[1]):
        max3indices = np.argpartition(avg_refl_coeffs[i,j,:],-3)[-3:]
        R_maxes = (avg_refl_coeffs[i,j,max3indices[0]],avg_refl_coeffs[i,j,max3indices[1]],avg_refl_coeffs[i,j,max3indices[2]])
        lambda_maxes = (avg_wvls[max3indices[0]],avg_wvls[max3indices[1]],avg_wvls[max3indices[2]])
        interp = af.parabolic_max(lambda_maxes,R_maxes)
        interp_peak_wvls[j,i] = interp[0].item()
        interp_max_Rs[j,i] = interp[1].item()
        
#lazy norm
interp_max_Rs = np.divide(np.abs(interp_max_Rs),np.amax(np.abs(interp_max_Rs)))
integrated_Rs = np.divide(np.abs(integrated_Rs),np.amax(np.abs(integrated_Rs)))

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.plot_surface(aa,bb,integrated_Rs,cmap="inferno",edgecolor="None")
ax.set_xlabel("α")
ax.set_ylabel("β")
ax.set_zlabel("integrated refl power")
if export_graphs:
    plt.savefig(f"{fig_output_directory}/integral-surface.pdf",dpi=500)
plt.show()

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.plot_surface(aa,bb,interp_max_Rs,cmap="viridis",edgecolor="None")
ax.set_xlabel("α")
ax.set_ylabel("β")
ax.set_zlabel("parabolically interpolated max refl power")
if export_graphs:
    plt.savefig(f"{fig_output_directory}/interp-surface.pdf",dpi=500)
plt.show()

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.plot_surface(aa,bb,interp_peak_wvls,cmap="hot",edgecolor="None")
ax.set_xlabel("α")
ax.set_ylabel("β")
ax.set_zlabel("interpolated intense wavelength (μm)")
if export_graphs:
    plt.savefig(f"{fig_output_directory}/wvl-surface.pdf",dpi=500)
plt.show()

# %%
