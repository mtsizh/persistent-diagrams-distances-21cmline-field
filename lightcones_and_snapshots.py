#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 17:09:28 2026
@author: maks
"""
import numpy as np 
import matplotlib.pyplot as plt
if "EoR" in plt.colormaps():
    plt.colormaps.unregister("EoR")
import py21cmfast as p21c 
import tempfile
p21c.config.USE_CACHE = False  
p21c.config.KEEP_DEFAULT_OUTPUT = False 


astro_params = {"F_STAR10": -1.3, # star formation efficiency (atomic cooling galaxies) for pivot mass 1e10 Msun (log10)
                "ALPHA_STAR": 0.5, # slope of the dependency of star formation efficiency on the host halo mass 
                "F_ESC10": -1.0, # escape fraction of Lyman photons into the IGM for pivot mass 1e10 Msun (log10)
                "ALPHA_ESC": -0.5, # slope of the dependency of escape fraction on the host halo mass 
                "L_X": 40.5, # X-ray luminosity (log10)
                "M_TURN": 8.7,
                "t_STAR": 0.5,
}
astro_options = {"USE_MINI_HALOS": False, "USE_UPPER_STELLAR_TURNOVER": False, 
                 "USE_EXP_FILTER": False, "USE_TS_FLUCT":True } 

nodes_array = np.logspace(0.7, 1.6, 25)
seeds = [111,222,333,444,555,666,777]
Om = [0.287,0.301,0.315,0.329,0.343]
#pwindx = [0.96,0.965,0.97] 
#sigma8 = [0.787,0.799,0.811,0.823,0.835]
#h0 = [0.67,0.68,0.69,0.70,0.71,0.72,0.73]

z_min = 7.0
z_max = 30.0

for i in seeds:
    for j in Om:
        with tempfile.TemporaryDirectory() as tmp_dir:
            cache1 = p21c.OutputCache(tmp_dir)
            cosmo_params1 = p21c.CosmoParams(
                SIGMA_8 = 0.811,
                hlittle=0.67, 
                OMb=0.0458, 
                OMm=j, 
                POWER_INDEX=0.965
            )
            inputs1 = p21c.InputParameters(random_seed=i, cosmo_params=cosmo_params1,
                                              node_redshifts=nodes_array.tolist(),
                                              astro_params=astro_params, astro_options=astro_options
            ).evolve_input_structs(BOX_LEN=300, DIM=300, HII_DIM=200)
            min_redshift=z_min
            max_redshift=z_max
            lcn1 = p21c.RectilinearLightconer.between_redshifts(
                min_redshift=min_redshift,
                max_redshift=max_redshift,
                resolution=inputs1.simulation_options.cell_size,
                cosmo=inputs1.cosmo_params.cosmo)
            lightcone1 = p21c.run_lightcone(
                lightconer=lcn1,
                inputs=inputs1,
                cache=cache1,
                regenerate=True)
            file_sim = str(i) +'_'+str(j) + ".h5"
            lightcone1.save(f"{i}_{j}.h5")
            print(f"Simulation saved in: {file_sim}") 
            '''
            lightcone_redshifts = lightcone1.lightcone_redshifts # numpy array
            idx_start = np.argmin(np.abs(lightcone_redshifts - z))
            map_at_z = lightcone1.lightcones['brightness_temp'][:, :, idx_start]
            file_snap = str(i) +'_'+str(j)
            np.save(file_snap, map_at_z)
            print(f"Snapshot saved: {file_snap}")
            '''
        
