'''
Python script to simulate the UHECR dataset for verification purposes.
This is converted from the notebook fit_to_simulation/run_simulation.ipynb.

This is made so that this can be run on command line as a bash script.
'''
import os
import numpy as np
import argparse

from fancy import Data, Model, Analysis
from fancy.interfaces.stan import get_simulation_input

# paths to important files
path_to_this_file = os.path.abspath(os.path.dirname(__file__))
stan_path = os.path.join(path_to_this_file, "stan")
source_file = os.path.join(path_to_this_file, "data", "sourcedata.h5")
table_path = os.path.join(path_to_this_file, "tables")
output_path = os.path.join(path_to_this_file, "output")

# make output path if it doesnt exist
if not os.path.exists(output_path):
    os.mkdir(output_path)

parser = argparse.ArgumentParser(description='Simulate UHECR dataset.')
parser.add_argument(
    '--source',
    dest='source_type',
    action='store',
    default="SBG_23",
    type=str,
    help='The source catalogue used (from SBG_23, 2FHL_250Mpc, swift_BAT_213)',
    choices=["SBG_23", "2FHL_250Mpc", "swift_BAT_213"])
parser.add_argument(
    '--detector',
    dest='detector_type',
    action='store',
    default="TA2015",
    type=str,
    help='The type of detector config (from TA2015, auger2014, auger2010)',
    choices=["TA2015", "auger2014", "auger2010"])
parser.add_argument(
    '--model',
    dest='model_type',
    action='store',
    default="joint",
    type=str,
    help='The stan model considered for simulation (from joint, joint_gmf)',
    choices=["joint", "joint_gmf"])
parser.add_argument(
    '--ptype',
    dest='ptype',
    action='store',
    default="p",
    type=str,
    help=
    "Type of particle used for back propagation (only used with joint_gmf model)."
)
parser.add_argument(
    "--seed",
    dest="seed",
    action="store",
    default=19990308,
    type=int,
    help="Random seed used for fixed parameter sampling in Stan.")
parser.add_argument("--sim_inputs",
                    dest="sim_inputs",
                    default=[0.5, 20, 3.0],
                    help="Simulation inputs: f, B, alpha",
                    nargs=3)

parser.add_argument("--tight_B",
                    dest="tight_B",
                    action="store_true",
                    default=False,
                    help="Use a tighter B-field of 1nG instead.")


def get_detectorimports(detector_type):
    '''Get variables imported by (detector_name).py'''
    if detector_type == "TA2015":
        from fancy.detector.TA2015 import detector_properties, alpha_T, M, Eth
    elif detector_type == "auger2014":
        from fancy.detector.auger2014 import detector_properties, alpha_T, M, Eth
    elif detector_type == "auger2010":
        from fancy.detector.auger2010 import detector_properties, alpha_T, M, Eth
    else:
        raise Exception("Undefined detector type!")

    return detector_properties, alpha_T, M, Eth


def get_sim_inputs(sim_inputs, D_src, M, alpha_T):
    '''Return parameters used for Model.input given simulation inputs.'''
    f, B, alpha = sim_inputs

    FT_PAO = 0.3601  # total flux using auger2014 data
    Nsim_expected = FT_PAO / (M / alpha_T)
    Nsim = int(np.round(Nsim_expected))

    # check value for Nsim
    print("Simulated events: {0}".format(Nsim))

    # L in yr^-1, F in km^-2 yr^-1
    L, F0 = get_simulation_input(Nsim, float(f), D_src, M, alpha_T)

    # check luminosity and isotropic flux values
    # L ~ O(10^39), F0 ~ 0.18
    # same luminosity so only need to check one value
    print("Simulated Luminosity: {0:.3e}".format(L[0]))
    print("Simulated isotropic flux: {0:.3f}".format(F0))

    return float(B), L, F0, float(alpha)


if __name__ == "__main__":

    args = parser.parse_args()

    # get filenames
    table_file = os.path.join(
        table_path, 'tables_{0}_{1}.h5'.format(args.source_type,
                                               args.detector_type))

    tightB_label = "tightB" if args.tight_B else "notightB"
    sim_output_file = os.path.join(
        output_path,
        "{0}_sim_{1}_{2}_{3}_{4}_{5}.h5".format(args.model_type,
                                                args.source_type,
                                                args.detector_type, args.seed,
                                                args.ptype, tightB_label))
    # sim_output_file = os.path.join(
    #     output_path,
    #     "{0}_sim_{1}_{2}_{3}_{4}.h5".format(args.model_type, args.source_type,
    #                                         args.detector_type, args.seed,
    #                                         args.ptype))

    # get things related to detector
    detector_properties, alpha_T, M, Eth = get_detectorimports(
        args.detector_type)

    # Create the Data() object
    data = Data()
    data.add_source(source_file, args.source_type)
    data.add_detector(detector_properties)

    # get the simulation inputs
    D_src = data.source.distance
    B, L, F0, alpha = get_sim_inputs(args.sim_inputs, D_src, M, alpha_T)
    Eth_sim = 20  # in EeV, set globally for now

    # for tighter B-field, use 1nG instead.
    if args.tight_B:
        B = 1.

    # create model, compile it, and add simulation inputs to it
    sim_name = os.path.join(stan_path,
                            '{0}_model_sim.stan'.format(args.model_type))

    simulation = Model(sim_filename=sim_name, include_paths=stan_path)
    simulation.compile(reset=False)

    simulation.input(B=B, L=L, F0=F0, alpha=alpha, Eth=Eth, ptype=args.ptype)

    # create Analysis object and perform the simulation (sampling)
    summary = b'Simulation of UHECR dataset'
    sim_analysis = Analysis(data,
                            simulation,
                            analysis_type=args.model_type,
                            filename=sim_output_file,
                            summary=summary)

    # building exposure and energy tables before sampling in stan
    sim_analysis.build_tables(sim_only=True)
    # simulate
    sim_analysis.simulate(seed=args.seed, Eth_sim=Eth_sim)
    # save file
    sim_analysis.save()

    # print resulting UHECR observed after propagation and Elosses
    print("Observed simulated UHECRs: {0}\n".format(
        len(sim_analysis.arrival_direction.unit_vector)))
