'''
Python script to simulate the UHECR dataset for verification purposes.
This is converted from the notebook fit_to_simulation/run_simulation.ipynb.

This is made so that this can be run on command line as a bash script.
'''
import os

from fancy import Data, Model, Analysis
from fancy.interfaces.stan import get_simulation_input

# paths to important files
path_to_this_file = os.path.abspath(os.path.dirname(__file__))
stan_path = os.path.join(path_to_this_file, "stan")
source_file = os.path.join(path_to_this_file, "data", "sourcedata.h5")
table_path = os.path.join(path_to_this_file, "tables")
output_path = os.path.join(path_to_this_file, "outputs")
output_sim_path = os.path.join(output_path, "sim")

# make output path if it doesnt exist
if not os.path.exists(output_path):
    os.mkdir(output_path)

if not os.path.exists(output_sim_path):
    os.mkdir(output_sim_path)


class Simulation():
    def __init__(self, sim_args, verbose=False):
        '''Perform simulation for some configuration.'''
        self.source = sim_args["source"]
        self.detector = sim_args["detector"]
        self.ptype = sim_args["ptype"]
        self.model = sim_args["model"]

        self.seed = sim_args["seed"]
        self.end_label = sim_args["end_label"]

        self.verbose = verbose

        self._init_files()

    def _init_files(self):
        '''Initialize filenames for output files / table files'''
        self.table_file = os.path.join(
            table_path, 'tables_{0}_{1}.h5'.format(self.source, self.detector))

        if self.end_label is not None:
            sim_output_file = "{0}_{1}_{2}_{3}_{4}.h5".format(
                self.source, self.detector, self.seed, self.ptype,
                self.end_label)
        else:
            sim_output_file = "{0}_{1}_{2}_{3}.h5".format(
                self.source, self.detector, self.seed, self.ptype)

        # make required directories
        if not os.path.exists(os.path.join(output_sim_path, self.model)):
            os.mkdir(os.path.join(output_sim_path, self.model))

        if not os.path.exists(
                os.path.join(output_sim_path, self.model, "simulation")):
            os.mkdir(os.path.join(output_sim_path, self.model, "simulation"))

        self.sim_output_file = os.path.join(output_sim_path, self.model,
                                            "simulation", sim_output_file)

        if self.verbose:
            print("Table output file: ", self.table_file)
            print("Simulation Output File: ", self.sim_output_file)

    def _get_detectorimports(self):
        '''Get variables imported by (detector_name).py'''
        if self.detector == "TA2015":
            from fancy.detector.TA2015 import detector_properties, alpha_T, M, Eth
        elif self.detector == "auger2014":
            from fancy.detector.auger2014 import detector_properties, alpha_T, M, Eth
        elif self.detector == "auger2010":
            from fancy.detector.auger2010 import detector_properties, alpha_T, M, Eth
        else:
            raise Exception("Undefined detector type!")

        return detector_properties, alpha_T, M, Eth

    def set_simulation(self, f=0.5, B=20, alpha=3, Nsim=2500):
        '''Setup the simulation'''
        print("Setting up the simulation...")
        detector_properties, alpha_T, M, Eth = self._get_detectorimports()

        # Create the Data() object
        data = Data()
        data.add_source(source_file, self.source)
        data.add_detector(detector_properties)

        # setup the simulation inputs
        print("Simulated events: {0}".format(Nsim))

        # L in yr^-1, F in km^-2 yr^-1
        L, F0 = get_simulation_input(Nsim, f, data.source.distance, M, alpha_T)

        # check luminosity and isotropic flux values
        # L ~ O(10^39), F0 ~ 0.18
        # same luminosity so only need to check one value
        print("Simulated Luminosity: {0:.3e}".format(L[0]))
        print("Simulated isotropic flux: {0:.3f}".format(F0))

        if self.verbose:
            print(
                "Inputs: Nsim = {0:d}, f = {1:.1f}, B = {2:.2f} nG, alpha = {3:d}, \n L = {4:.3e} 10^-38 yr^-1, F0 = {5:.3f} km^-2 yr^-1"
                .format(Nsim, f, B, alpha, L[0] * 1e-38, F0))

        # create model, compile it, and add simulation inputs to it
        sim_name = os.path.join(stan_path,
                                '{0}_model_sim.stan'.format(self.model))

        simulation = Model(sim_filename=sim_name, include_paths=stan_path)
        simulation.compile(reset=False)

        simulation.input(B=B,
                         L=L,
                         F0=F0,
                         alpha=alpha,
                         Eth=Eth,
                         ptype=self.ptype)

        if self.verbose:
            print("Model used for simulation: ", sim_name)

        summary = b'Simulation of UHECR dataset'
        self.sim_analysis = Analysis(data,
                                     simulation,
                                     analysis_type=self.model,
                                     filename=self.sim_output_file,
                                     summary=summary)

        # building exposure and energy tables before sampling in stan
        self.sim_analysis.build_tables(sim_only=True)

    def simulate(self, Eth_sim=20):
        '''Perform simulation'''
        if self.verbose:
            print("Performing Simulation with seed: ", self.seed)
        self.sim_analysis.simulate(seed=self.seed, Eth_sim=Eth_sim)

        # print resulting UHECR observed after propagation and Elosses
        print("Observed simulated UHECRs: {0}\n".format(
            len(self.sim_analysis.arrival_direction.unit_vector)))

    def save(self):
        '''Save to simulation output file'''
        self.sim_analysis.save()
