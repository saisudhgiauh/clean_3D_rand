# gen_data_standard.py
# standard sequential version
# - Same k-point convention (scales with f_min, frequency-dependent angle mapping)
# - Same source placement and components
# - Same PML and eps_averaging settings (taken from config)
# - Same near-to-far setup
# - Saves wvls, inc_angles, incident farfield dict (inc), reflected farfield dict (refl)
# - Uses polarisation-aware stop condition via ff.pol_comp(...)

import os

import numpy as np
import meep as mp


# Local project imports
import geometry as gg
import libs.field_functions as ff
from config import *  # expects: cell_x, cell_y, cell_z, pml_layers, resolution, n_freq,
                      #          f_mu, f_sigma, f_min, n2f_width, n2f_height, n2f_res,
                      #          source_y, incident_variation, n_ridges, max_ridge_angle,
                      #          randomise_ridge_angles, n_periods, eps_averaging,
                      #          air, data_output_directory, bloch_angles,
                      #          polarisation_setting






import time





def _build_geometry():
     """Reconstruct grating/ridge geometry exactly as in the sequential version."""
     if randomise_ridge_angles:
        etas = np.load("rands.npy")
     else:
         etas = [max_ridge_angle for _ in range(n_ridges)]

     actual_geometry = []
     for i in range(len(etas)):
        position = (2 * i + 1 - n_ridges) / (2 * n_ridges) * cell_x
        actual_geometry += gg.get_ridge(position, etas[i])
     return actual_geometry



def _make_sources(polarisation: str):
    sources = []
    if polarisation == "P":
        sources.append(
            mp.Source(
                mp.GaussianSource(f_mu, fwidth=f_sigma),
                component=mp.Ex,
                center=mp.Vector3(0, source_y, 0),
                size=mp.Vector3(cell_x-2*pml_width, 0, cell_z-2*pml_width),
            )
        )
    elif polarisation == "S":
        sources.append(
            mp.Source(
                mp.GaussianSource(f_mu, fwidth=f_sigma),
                component=mp.Ez,
                center=mp.Vector3(0, source_y, 0),
                size=mp.Vector3(cell_x-2*pml_width, 0, cell_z-2*pml_width),
            )
        )
    else:  # mixed
        sources.append(
            mp.Source(
                mp.GaussianSource(f_mu, fwidth=f_sigma),
                component=mp.Ex,
                center=mp.Vector3(0, source_y, 0),
                size=mp.Vector3(cell_x, 0, cell_z),
            )
        )
        sources.append(
            mp.Source(
                mp.GaussianSource(f_mu, fwidth=f_sigma),
                component=mp.Ez,
                center=mp.Vector3(0, source_y, 0),
                size=mp.Vector3(cell_x, 0, cell_z),
            )
        )
    return sources


def _k_point_from_bloch_angle(bloch_angle_deg: float) -> mp.Vector3:
    """Match sequential script: use f_min scaling and incident_variation switch."""
    ba_rad = np.deg2rad(bloch_angle_deg)
    if incident_variation == "alpha":
        return mp.Vector3(np.sin(ba_rad), -1 * np.cos(ba_rad), 0).scale(f_min)
    else:
        return mp.Vector3(0, -1 * np.cos(ba_rad), np.sin(ba_rad)).scale(f_min)



def run_angle(angle_deg: float, polarisation: str) -> str:
   





    # --- Boundary layers & near-to-far surface (create locally to avoid pickling) ---
    pml_layers = [mp.PML(pml_width)]
   


    """Single-angle worker that reproduces sequential physics & outputs."""
    # Skip if already completed (idempotent)
    done_marker = os.path.join(data_output_directory, f"done-{angle_deg}.txt")
    if os.path.exists(done_marker):
        return f"Skipping {angle_deg} — already done"

    # Geometry and detectors
    actual_geometry = _build_geometry()

    n2f_surface = mp.Near2FarRegion(
        center=mp.Vector3(0, (cell_y / 2) - (2 * pml_width), 0),
        size=(cell_x , 0, cell_z )
    )

    # Sources & k-point
    sources = _make_sources(polarisation)
    wavevector = _k_point_from_bloch_angle(angle_deg)

    # --- Normalisation run ---
    normalisation_geometry = [
        mp.Block(
            mp.Vector3(cell_x, cell_y, cell_z),
            center=mp.Vector3(),
            material=air,  # match sequential behavior (not strict vacuum)
        )
    ]

    norm_sim = mp.Simulation(
        cell_size=mp.Vector3(cell_x, cell_y, cell_z),
        boundary_layers=pml_layers,
        geometry=normalisation_geometry,
        sources=sources,
        resolution=resolution,
        k_point=wavevector,
        eps_averaging=eps_averaging,
    )

    inc_n2f = norm_sim.add_near2far(f_mu, f_sigma, n_freq, n2f_surface)

    norm_sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            10,
            ff.pol_comp(polarisation),
            mp.Vector3(0, (cell_y / 2) - (2 * pml_width), 0),
            1e-6,
        )
    )

    

    # Save incident farfield at the reference point (for parity with sequential script)
    farfield_source = norm_sim.get_farfields(
        inc_n2f,
        n2f_res,
        center=mp.Vector3(0, n2f_height, 0),
        size=mp.Vector3(),
    )

    # Near-to-far subtraction payload
    farfield_source_subtractme = norm_sim.get_near2far_data(inc_n2f)

    # --- Actual run with geometry ---
    # Important to reset Meep before the new Simulation
    norm_sim.reset_meep()

    sim = mp.Simulation(
        cell_size=mp.Vector3(cell_x, cell_y, cell_z),
        boundary_layers=pml_layers,
        geometry=actual_geometry,
        sources=sources,
        resolution=resolution,
        k_point=wavevector,
        eps_averaging=eps_averaging,
    )

    refl_n2f = sim.add_near2far(
        f_mu, f_sigma, n_freq, n2f_surface, nperiods=n_periods
    )
    sim.load_minus_near2far_data(refl_n2f, farfield_source_subtractme)

    sim.run(
        until_after_sources=mp.stop_when_fields_decayed(
            10,
            ff.pol_comp(polarisation),
            mp.Vector3(0, (cell_y / 2) - (2 * pml_width), 0),
            1e-6,
        )
    )

    # Collect reflected farfield on grid (same as sequential)
    farfield_ridge = sim.get_farfields(
        refl_n2f,
        n2f_res,
        center=mp.Vector3(0, n2f_height, 0),
        size=mp.Vector3(n2f_width, 0, n2f_width),
    )

    # Frequencies / wavelengths and incident angles per frequency
    freqs = mp.get_near2far_freqs(refl_n2f)
    wavelengths = np.reciprocal(freqs)

    inc_angles = np.zeros(n_freq)
    if incident_variation == "alpha":
        for i in range(len(inc_angles)):
            inc_angles[i] = np.rad2deg(np.arcsin(wavevector.x / freqs[i]))
    else:
        for i in range(len(inc_angles)):
            inc_angles[i] = np.rad2deg(np.arcsin(wavevector.z / freqs[i]))

    # --- Persist exactly the same outputs as the sequential script ---
    os.makedirs(data_output_directory, exist_ok=True)
    np.save(f"{data_output_directory}/wvls-{angle_deg}", wavelengths)
    np.save(f"{data_output_directory}/inc_angles-{angle_deg}", inc_angles)
    np.save(f"{data_output_directory}/inc-{angle_deg}", farfield_source)
    np.save(f"{data_output_directory}/refl-{angle_deg}", farfield_ridge)

    # Marker file
    open(done_marker, "w").close()

    return f"Done ({angle_deg})"


for ba in bloch_angles:
    run_angle(ba, polarisation_setting)



