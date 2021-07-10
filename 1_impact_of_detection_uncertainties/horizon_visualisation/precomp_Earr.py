#!/usr/bin/env python
from mpi4py import MPI
import numpy as np
import h5py
import time

from fancy.propagation.energy_loss import get_arrival_energy

"""
Script to launch parallel calculation of UHECR propagation losses 
using the continuous loss approximation.

Store the resulting Earr in output file.
"""

COMM = MPI.COMM_WORLD

if COMM.rank == 0:
    
    start_time = time.time()

    # Define output
    output_file = 'output/my_precomp_Earr.h5'

    # Set parameters
    Eth = 50
    Emax = 1.0e5
    N = 100
    E_grid = np.logspace(np.log(Eth), np.log(Emax), N, base = np.e)
    D_grid = np.linspace(0, 500, N)
    
    Es = np.array_split(E_grid, COMM.size)

    # Store parameters
    with h5py.File(output_file, 'w') as f:
        f.create_dataset('Eth', data = Eth)
        f.create_dataset('Emax', data = Emax)
        f.create_dataset('E_grid', data = E_grid)
        f.create_dataset('D_grid', data = D_grid)
        f.create_dataset('N', data = N)

    D = D_grid
    
else:
    
    Es = None
    D = None

Es = COMM.scatter(Es, root=0)
D = COMM.bcast(D, root=0)

Earr = np.zeros((len(Es), len(D)))

# Compute the arrival energies
for i, E in enumerate(Es):
    
    for j, d in enumerate(D):
        
        Earr[i][j] = get_arrival_energy(E, d)

    print(i, E, 'completed')
    
results = MPI.COMM_WORLD.gather(Earr, root=0)

if COMM.rank == 0:
    results = np.vstack(results)
    print('time:', time.time() - start_time)

    with h5py.File(output_file, 'r+') as f:
        f.create_dataset('Earr_grid', data = results)
        
    
