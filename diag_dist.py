#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 11:09:00 2026
@author: maks
"""
import numpy as np 
import matplotlib.pyplot as plt
if "EoR" in plt.colormaps():
    plt.colormaps.unregister("EoR")
import py21cmfast as p21c 
import gudhi as gd
import gudhi.wasserstein
import os
import pandas as pd
from skbio.stats.distance import DistanceMatrix, permanova

zlist=[7,7.5,8,8.5,9,9.5,10,10.5,11,10.5,11,11.5,12,12.5,13,13.5,14]
labelslist = []
diagl = []

def load_lightcone(seed, om, data_dir="."):
    fname = os.path.join(data_dir, f"{seed}_{om}.h5")
    if not os.path.exists(fname):
        raise FileNotFoundError(fname)
    return p21c.LightCone.from_file(fname)

seeds = [111,222,333,444,555,666,777]
cosmology = [0.67,0.72]
data_dir = "300h0/data"   
lightcones = {}
for o in cosmology:
    for s in seeds:
        try:
            lightcones[(s, o)] = load_lightcone(s, o, data_dir)
        except FileNotFoundError as e:
            print(f"[warning] missing file, skipping: {e}")
if not lightcones:
    raise SystemExit("No lightcone files found -- check `data_dir`, `seeds`, `cosmology`.")

for (o, s), lc in lightcones.items():
    loaded_lightcone = lc
    lightcone_redshifts = loaded_lightcone.lightcone_redshifts # numpy array
    labelslist.append(str(s)+'_'+str(o))
    diaglistforz = []
    for z in zlist:
            idx_start = np.argmin(np.abs(lightcone_redshifts - z))
            redshifts = np.array(loaded_lightcone.node_redshifts)
            map_at_z = loaded_lightcone.lightcones['brightness_temp'][:, :, idx_start]
            map_at_z = map_at_z - np.mean(map_at_z)
            cc = gd.CubicalComplex(dimensions=map_at_z.shape, top_dimensional_cells=map_at_z.flatten())    
            diag = cc.persistence()
            diag = [item[1] for item in diag if item[0] == 1]
            diag = np.array(diag)
            diaglistforz.append(diag)
    diagl.append(diaglistforz)
print('All file successfully read, diagrams computed')
    
plist = []
maxlist = []
for k, z in enumerate(zlist):
        n = len(labelslist)
        dist_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                dist = gudhi.wasserstein.wasserstein_distance(diagl[i][k], diagl[j][k], order=2.0)
                dist_matrix[i, j] = dist
                dist_matrix[j, i] = dist 
                print('Wasserstein distance for z=',z, ';i=', i,'j=',j,'computed')
                    
        fig, ax = plt.subplots(figsize=(8, 8))
        im = ax.imshow(dist_matrix, cmap='viridis')
        cbar = ax.figure.colorbar(im, ax=ax)
        cbar.ax.set_ylabel("Wasserstein Distance", rotation=-90, va="bottom")
        ax.tick_params(axis='x', labelrotation=45)
        ax.set_xticks(np.arange(n))
        ax.set_yticks(np.arange(n))
        ax.set_xticklabels(labelslist)
        ax.set_yticklabels(labelslist)
        
        ids = [f'sim_{i}' for i in range(len(cosmology)*len(seeds))]
        dm = DistanceMatrix(dist_matrix, ids)
        maxlist.append(np.max(dist_matrix))
        metadata = pd.DataFrame({
            'cosmo_model': ['Group_A']*7 + ['Group_B']*7  
#            ['Group_C']*5 + ['Group_D']*5 + 
#            ['Group_E']*5 + ['Group_F']*5
        }, index=ids)
        results = permanova(dm, metadata, column='cosmo_model', permutations=1999)
        p_value = round(results['p-value'],7)
        plist.append(p_value)
        print('z=',z,' p=', p_value)
        f_stat = results['test statistic']

        titlstr="Wasserstein Distances, z="+str(z) + ", p-value="+str(p_value)
        ax.set_title(titlstr)
        fig.tight_layout()
        plt.savefig('dist_mat_h0'+str(z)+'.png', dpi=600)
        plt.show()
with open('plist.txt', 'w') as f:
    for p in plist:
        f.write(f"{p}\n")
with open('maxkist.txt', 'w') as f:
    for p in maxlist:
        f.write(f"{p}\n")
