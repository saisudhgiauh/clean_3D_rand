import meep as mp
import numpy as np

def pol_comp(polarisation_setting):

    '''
    inputs:
        polarisation_setting (string, "P", "S" or "mixed") -- the polarisation config variable
    outputs:
        pol_comp (int, alias for meep field component) -- the MEEP field component corresponding to the polarisation config variable.
            If not a single polarisation, mp.Ex is returned arbitrarily. Used for terminating MEEP sims after sufficient field decay.
    '''

    if polarisation_setting == "P":
        return mp.Ex
    elif polarisation_setting == "S":
        return mp.Ez
    else:
        return mp.Ex
    
def pol_string(polarisation_setting):

    '''
    inputs: 
        polarisation_setting (string, "P", "S", or "mixed") -- the polarisation config variable
    outputs:
        pol_string (string) -- the field component (in string form) corresponding to the polarisation config variable. If not a single polarisation,
            "Ex" is returned arbitrarily. Used for extracting field data from n2f dictionaries with keys like "Ex", "Ey", "Ez".
    '''

    if polarisation_setting == "P":
        return "Ex"
    elif polarisation_setting == "S":
        return "Ez"
    else:
        return "Ex"
    
def h5_pol_string(polarisation_setting):

    '''
    inputs: 
        polarisation_setting (string, "P", "S" or "mixed") -- the polarisation config variable
    outputs:    
        h5_pol_string (string) -- the REAL field component (in string form) corresponding to the polarisation config variable. If not a single polarisation,
            "ex.r" is returned arbitrarily. Used for extracting field data from datasets within .h5 files. It is recommended to use MEEP to write to numpy
            arrays and handle these via np.save and np.load instead of using the h5 format, which sometimes has problems when used with MPI parallelisation.
    '''

    if polarisation_setting == "P":
        return "ex.r"
    elif polarisation_setting == "S":
        return "ez.r"
    else:
        return "ex.r"
    
'''an unused function that calculates the farfield radial poynting vector given farfield parameters and an angle, which was superceded by direct generation
of angular data via n2f grids.

def poynting(angle,y_distance,sim,n2f_surface,n_freq):
    E = np.zeros((n_freq,3),dtype=np.complex128)
    H = np.zeros((n_freq,3),dtype=np.complex128)

    detector_coords = mp.Vector3(y_distance * np.tan(np.deg2rad(angle)),y_distance,0)
    far_field = np.array(sim.get_farfield(n2f_surface,detector_coords))
    split = np.array_split(far_field,2*n_freq)
    for i in range(2*n_freq):
        if i%2 == 0: #i even, electric field component
            E[int(i/2),:] = np.conj(split[i])
        else:   #i odd, magnetic field component
            H[int((i-1)/2),:] = split[i]

    Px = np.real(E[:,1]*H[:,2]-E[:,2]*H[:,1]) # Ey*Hz-Ez*Hy
    Py = np.real(E[:,2]*H[:,0]-E[:,0]*H[:,2]) # Ez*Hx-Ex*Hz
    Pr = np.sqrt(np.square(Px)+np.square(Py))

    return(
        Pr,
        E,
        H
    )
'''
    
    
