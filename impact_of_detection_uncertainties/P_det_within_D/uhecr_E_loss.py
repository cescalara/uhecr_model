#!/usr/bin/env python

import numpy as np
from mpi4py import MPI
import h5py
import time
import sys

import stan_utility


"""
Script to launch parallel calculation of UHECR energy losses for the calculation
of P(E > E_th | D) and P(\hat{E} > E_th | D).

Variable names:
E -> Earr
\hat{E} -> Edet
E_th -> Eth
"""


# Simulation parameters
alpha = 2.0
Eth_sim = 20
Eth = 52.0

# Number of UHECR per simulation
Ncr = 100 # 500, 1000

# Number of simulations at each D
Ntrials = 4 # 48*10 

# Range of distances simulated
Ds = np.linspace(0, 500, 10) # np.linspace(0, 500, 100) 

# Define Stan simulation model
sim_filename = 'uhecr_E_loss.stan'

# Define output HDF5 file 
output_file = 'simulation_output/my_uhecr_E_loss_output.hdf5'


def run_stan_sim(N, Eth_sim, alpha, D, Eth, sim_filename):
    """
    Run the Stan simulation for N events above Eth_sim from distance D
    and return the fraction above Eth.
    
    :param N: Number of UHECRs to simulate.
    :param Eth_sim: The minimum energy of UHECR generated in the simulation.
    :param alpha: The spectral index of the sourc UHECR (power law spectrum).
    :param D: The distance at which a shell of sources is to be placed.
    :param Eth: The threshold energy of a UHECR sample.
    :param sim_filename: The filename of the Stan simulation code.
    
    :return: The detection probability for Edet > Eth and Earr > Eth.
    """

    # Run the simulation.
    sim_input = {'N' : N, 'alpha' : alpha, 'Eth_sim' : Eth_sim, 'D' : D}
    sim = stan_utility.compile_model(filename = sim_filename, model_name = 'uhecr_E_loss')
    sim_output = sim.sampling(data = sim_input, iter = 1, chains = 1, 
                              algorithm = "Fixed_param")

    # Extract the output.
    E = sim_output.extract(['E'])['E'][0]
    Earr = sim_output.extract(['Earr'])['Earr'][0]
    Edet = sim_output.extract(['Edet'])['Edet'][0]

    # Count number above threshold
    N_arr_gt_Eth = np.shape(np.where(Earr > Eth))[1]
    N_det_gt_Eth = np.shape(np.where(Edet > Eth))[1]

    return (N_arr_gt_Eth / N), (N_det_gt_Eth / N)


COMM = MPI.COMM_WORLD


if COMM.rank == 0:

    # Initialise output file and compile stan code to cache 

    with h5py.File(output_file, 'w') as f:
        f.create_dataset('Parr', (len(Ds),), 'f') 
        f.create_dataset('Pdet', (len(Ds),), 'f')
        f.create_dataset('D', (len(Ds),), 'f')

    sim = stan_utility.compile_model(filename = sim_filename, model_name = 'uhecr_E_loss')

    done = False

    
# Start loop over Ds
for i, d in enumerate(Ds):
  
    if COMM.rank == 0:

        #print('D:', d)

        start_time = time.time()

        
        Eth = Eth
        Ncr = Ncr 
        alpha = alpha
        Eth_sim = Eth_sim
        D = d
        i = i
        sim_filename = sim_filename
        
        trial = range(Ntrials)
        trials = np.array_split(trial, COMM.size)
    
    else:

        trials = None
        N = None
        alpha = None
        Eth_sim = None
        Eth = None
        D = None
        sim_filename = None
        done = None

        
    # Distribute information
    trials = COMM.scatter(trials, root=0)
    Ncr = COMM.bcast(Ncr, root=0)
    alpha = COMM.bcast(alpha, root=0)
    Eth_sim = COMM.bcast(Eth_sim, root=0)
    Eth = COMM.bcast(Eth, root=0)
    D = COMM.bcast(D, root=0)
    sim_filename = COMM.bcast(sim_filename, root=0)
    done = COMM.bcast(done, root=0)
    
    Parr = np.zeros(len(trials))
    Pdet = np.zeros(len(trials))
    
    # Run parallel trials for each D value
    for j, t in enumerate(trials):

        if not done:
            
            try:

                Parr[j], Pdet[j] = run_stan_sim(Ncr, Eth_sim, alpha, D, Eth, sim_filename)

            except:

                Parr[j] = np.nan
                Pdet[j] = np.nan

        else:

            Parr[j] = 0.0
            Pdet[j] = 0.0
            
    # Gather results
    Parrs = MPI.COMM_WORLD.gather(Parr, root=0)
    Pdets = MPI.COMM_WORLD.gather(Pdet, root=0)

    if COMM.rank == 0:

        # Calculate mean
        Parr = np.nanmean([item for sublist in Parrs for item in sublist])
        Pdet = np.nanmean([item for sublist in Pdets for item in sublist])

        print('Parr:', Parr)
        print('Pdet:', Pdet)
        print('time:', time.time() - start_time)

        # Append results to file
        with h5py.File(output_file, 'r+') as f:
            f['Parr'][i] = Parr
            f['Pdet'][i] = Pdet
            f['D'][i] = D
            check_arr = f['Pdet'].value[i-5:i]

        # If Pdet has been 0 for a while, stop simulating.    
        if np.nanmean(check_arr) == 0.0:

            done = True
