import meep as mp
import numpy as np

#########################
###editable parameters###
#########################

###geometric parameters, all distances are in microns, all angles are in degrees.

basal_height =  0.395
branch_height = 0.06
branch_width = 0.28                 #in the case of tapering, this is the MAXIMUM width
trunk_width = 0.06
branch_sep = 0.14
tree_sep = 0.13
trunk_height = 1.46                 
n_branches = 14                     #may break if an odd number is used with offset_branches = False
start_taper_at = 10                 #branch number to start the taper at, measuring from the base (the one after this will be the first shortened branch). if no branch offset, must be an even number.
taper_level = 0.75                  #controls the aggressiveness of the taper, [0,1]. 1 means the top branch is nonexistent. 0 means no taper.
offset_branches = True             #whether or not the branches should be placed in an alternating offset pattern
round_edges = True                  #whether or not the ends of branches should be smoothed
round_centres = True                #whether or not the space between branches should be smoothed
round_top = True                    #whether or not the top of the tree (and overlap, if 3D) should be smoothed
point_spacing = 2                   #the distance along the ridge between overlapping layers
overlap = 0.5                       #the distance along the ridge that the overlap takes up
overlap_angle = 2                   #the angle of the overlapping part relative to the main ridge. For no overlap, set to something small but nonzero.
n_ridges = 4                        #the number of ridges that should be drawn in a single unitcell. Should always be 1 if randomise_ridge_angles is False.
randomise_ridge_angles = True      #whether or not the relative angles of the ridges should be randomised.
max_ridge_angle = 30                #the maximum angle from the normal that a ridge may be oriented at. If no randomisation, all ridges will be drawn at this angle.

###simulation parameters

pml_width = 0.1                     #thickness of the PML planes orthogonal to the y axis
refractive_index = 1.56 + 0.06j     #refractive index of the non-air parts of the unit cell
n_freq = 300                        #number of frequencies to sample with each DFT monitor.

f_min = 1/0.7                      #~min frequency, used to calculate stdev of input distribution (this will not be the true minimum)
f_max = 1/0.3                      #~max frequency, used to calculate stdev of input distribution (this will not be the true maximum)
f_mu = (f_min+f_max)/2                      #average of gaussian frequency input
resolution = 50                     #pixels per micron
n2f_height = 1e3                    #far field distance in the y direction (angular data calculated on this line(2D)/plane(3D))
n_periods = 6                      #number of periods simulated for the far field calculation (this is rather cpu intensive in 3D)
max_angle = 20                      #maximum angle away from the normal for the simulation to calculate the farfields.
n_angles = 61                       #number of angles in each direction to sample, set to an odd value so there is a true normal measurement
polarisation_setting = "mixed"      #set to "P" or "S" or "mixed". P polarisation has an electric field in the x direction, S polarisation has an electric field in the z direction
eps_averaging = False               #whether or not the dielectric should be eps averaged. sometimes takes a long time and I'm not sure why. Doesn't really affect the spectrum.
incident_variation = "beta"         #set to "alpha" or "beta". If multiple bloch angles are being used, this controls the direction in which the incident angle is varied between sims.
min_bloch_angle = 0                 #bloch wavevector angle of the first sim.
max_bloch_angle = 50                #bloch wavevector angle of the last sim.
bloch_angle_interval = 2            #change in bloch wavevector angle between sims

###output parameters

export_graphs = True               #whether or not any graphs produced should be exported (as pdf)
output_directory = "./visualiser"
data_output_directory = "./data"         #directory in which any data outputs (.npy files) should be placed. Must be pre-existing.
fig_output_directory = "./figures"       #directory in which any figures and spectra should be placed. Must be pre-existing.
data_input_directory = "./data"          #directory from which any data inputs (.npy files) should be accessed for plotting.

'''
inc_slice_min = 0
inc_slice_max = 20
inc_slice_interval = 5
'''

#############################
###non-editable parameters###
#############################

#stdev of source distribution
f_sigma = (f_max - f_min)
#air material
air = mp.Medium(epsilon = 1.0006)
#wing material. imaginary part of refractive index is handled by making the wing conductive. see MEEP docs for further info.
wing_mat = mp.Medium(epsilon = np.real(refractive_index)**2, D_conductivity = 2 * np.pi * (1/0.5) * 2* np.imag(refractive_index)/np.real(refractive_index))
#cell dimensions, obtained automatically from geometric parameters. currently the space at the top is automatically equal to the basal height.
cell_x = n_ridges * 2 * ((trunk_width / 2) + branch_width + (tree_sep / 2))
cell_y = 2 * basal_height + trunk_height 
cell_z = point_spacing
cell = mp.Vector3(cell_x,cell_y,cell_z)
#perfectly matched layers to truncate simulation in the Y direction beyond the extent of the cell
pml_layers = [mp.Absorber(pml_width,mp.Y)]
#y coordinate at which the planar source is drawn. currently 3/4 of the distance between the top of the trunk and the top of the cell.
source_y = (trunk_height / 2) + (3/4) * ((cell_y / 2) - (trunk_height / 2))
#x-z width of the square grid over which the farfields are calculated
n2f_width = 2*n2f_height * np.tan(np.deg2rad(max_angle))
#resolution of the near2far grid
n2f_res = n_angles / n2f_width
#array of bloch angles to perform separate simulations for
bloch_angles = np.arange(min_bloch_angle,max_bloch_angle + bloch_angle_interval,bloch_angle_interval)
