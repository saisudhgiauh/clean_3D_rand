#%%
import meep as mp
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0,os.getcwd())   #replace this with local directory
from config import *
import geometry as gg
import libs.field_functions as ff
import time

start_time = time.time()

actual_geometry = []
    
###construct geometry from config
    
if randomise_ridge_angles:
    etas = np.load("rands.npy")
else:
    etas = [max_ridge_angle for _ in range(n_ridges)]
for i in range(len(etas)):
    position = (2*i + 1 - n_ridges) / (2 * n_ridges) * cell_x
    actual_geometry += gg.get_ridge(position,etas[i])

###placing detectors and setting up normalisation geometry

n2f_surface = mp.Near2FarRegion(
    center=mp.Vector3(0,(cell_y / 2) - (2*pml_width),0),
    size=1.01 * mp.Vector3(cell_x,0,cell_z)    #transform plane surface
)

normalisation_geometry = [
        mp.Block(
            mp.Vector3(cell_x,cell_y,cell_z),    #block the same size as the cell, might be a better way of doing this but I can't find it
            center=mp.Vector3(),
            material=air #<-- this can't be a vacuum (epsilon=1) for some stupid reason
        )
    ]

###function to generate the spectrum for a single polarisation setting

def calc_single_spectrum(polarisation,bloch_angle):

    '''
    inputs:
        polarisation (string, "P", "S" or "mixed")
        bloch_angle (0 < float < 90) -- approx angle of propagation of minimum frequency component in simulation. The angle of incidence is
            frequency dependent - see pages 83-85 of 'Advances in FDTD Computational Electrodynamics : Photonics and Nanotechnology' for theory.
    outputs:
        wavelengths (1D array, length n_freq) -- the frequencies picked up by the DFT
        inc_angles (1D array, length n_freq) -- the incident propagation angles of each discrete frequency component in the input
        farfield_source (dictionary of 1D arrays, each length n_freq) -- incident electric field data at the central point of the n2f grid. Keys are "Ex", "Ey", "Ez" etc.
        farfield_ridge (dictionary of 3D arrays, each n_angles x n_angles x n_freq) -- reflected electric field data across the n2f grid.
    '''


    ba_rad = np.deg2rad(bloch_angle)
    wavevector = mp.Vector3(np.sin(ba_rad),-1*np.cos(ba_rad),0).scale(f_min) if incident_variation == "alpha" else mp.Vector3(0,-1*np.cos(ba_rad),np.sin(ba_rad)).scale(f_min)
    sources = []

    if polarisation == "P":     #P polarisation
        sources.append(
            mp.Source(mp.GaussianSource(f_mu,fwidth=f_sigma),
                component=mp.Ex,
                center=mp.Vector3(0,source_y,0), 
                size=mp.Vector3(cell_x,0,cell_z)
            )
        )
    elif polarisation == "S":  #S polarisation
        sources.append(
            mp.Source(mp.GaussianSource(f_mu,fwidth=f_sigma),
                component=mp.Ez,
                center=mp.Vector3(0,source_y,0), 
                size=mp.Vector3(cell_x,0,cell_z)
            )
        )
    elif polarisation == "mixed": #no polarisation
        sources.append(
            mp.Source(mp.GaussianSource(f_mu,fwidth=f_sigma),
                component=mp.Ex,
                center=mp.Vector3(0,source_y,0), 
                size=mp.Vector3(cell_x,0,cell_z)
            )
        ),
        sources.append(
            mp.Source(mp.GaussianSource(f_mu,fwidth=f_sigma),
                component=mp.Ez,
                center=mp.Vector3(0,source_y,0), 
                size=mp.Vector3(cell_x,0,cell_z)
            )
        )

    ###normalisation cell

    normalisation_sim = mp.Simulation(
        cell_size=mp.Vector3(cell_x,cell_y,cell_z),
        boundary_layers=pml_layers,
        geometry=normalisation_geometry,
        sources=sources,
        resolution=resolution,
        k_point=wavevector,
        eps_averaging=eps_averaging
    )

    ###first, we need the incident poynting vector flux with no messy geometry so we can 
    ###subtract it later and normalise against it
    print("Now performing normalisation run")

    inc_n2f = normalisation_sim.add_near2far(f_mu,f_sigma,n_freq,n2f_surface)
    normalisation_sim.run(
        until_after_sources=mp.stop_when_fields_decayed(10,ff.pol_comp(polarisation_setting),mp.Vector3(0,(cell_y / 2) - (2*pml_width),0),1e-6)
    )
    farfield_source = normalisation_sim.get_farfields(inc_n2f,n2f_res,center=mp.Vector3(0,n2f_height,0),size=mp.Vector3())
    print("Done with farfield norm")
    farfield_source_subtractme = normalisation_sim.get_near2far_data(inc_n2f)

    ###now for the interesting run with the christmas tree structure
    print("Now performing actual run")

    normalisation_sim.reset_meep()
    actual_sim = mp.Simulation(
        cell_size=mp.Vector3(cell_x,cell_y,cell_z),
        boundary_layers=pml_layers,
        geometry=actual_geometry,
        sources=sources,
        resolution=resolution,
        k_point=wavevector,
        eps_averaging=eps_averaging
    )
    refl_n2f = actual_sim.add_near2far(f_mu,f_sigma,n_freq,n2f_surface,nperiods=n_periods)
    actual_sim.load_minus_near2far_data(refl_n2f,farfield_source_subtractme)
    actual_sim.run(
        until_after_sources=mp.stop_when_fields_decayed(10,ff.pol_comp(polarisation_setting),mp.Vector3(0,(cell_y / 2) - (2*pml_width),0),1e-6)
    )
    farfield_ridge = actual_sim.get_farfields(refl_n2f,n2f_res,center=mp.Vector3(0,n2f_height,0),size=mp.Vector3(n2f_width,0,n2f_width))
    print("Done with farfield actual")

    ###things to spit out for plotting

    freqs = mp.get_near2far_freqs(refl_n2f)
    wavelengths = np.reciprocal(freqs)
    inc_angles = np.zeros(n_freq)
    if incident_variation == "alpha":
        for i in range(len(inc_angles)):
            inc_angles[i] = np.rad2deg(np.arcsin(wavevector.x/freqs[i]))
    else:
        for i in range(len(inc_angles)):
            inc_angles[i] = np.rad2deg(np.arcsin(wavevector.z/freqs[i]))
    return(
        wavelengths,
        inc_angles,
        farfield_source,
        farfield_ridge
    )

for ba in bloch_angles:
    output_prefix = f"{data_output_directory}/done-{ba}.txt"

    # Skip if already completed
    if os.path.exists(output_prefix):
        print(f"Skipping {ba} — already done")
        continue

    wvls,inc_angles,inc,refl = calc_single_spectrum(polarisation_setting,ba)    #perform separate simulations with each bloch angle.
    np.save(f"{data_output_directory}/wvls-{ba}",wvls)
    np.save(f"{data_output_directory}/inc_angles-{ba}",inc_angles)
    np.save(f"{data_output_directory}/inc-{ba}",inc)
    np.save(f"{data_output_directory}/refl-{ba}",refl)

    # Create a marker file so we know it's done
    open(output_prefix, "w").close()
    
    print(f"Done ({ba})")
# %%
