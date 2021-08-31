'''
Python script of the fitting process in the notebook run_simulation.ipynb.

This is made so that this can be run on command line as a bash script.
'''
import os
from fancy import Data, Model, Analysis

import argparse

# paths to important files
path_to_this_file = os.path.abspath(os.path.dirname(__file__))
stan_path = os.path.join(path_to_this_file, "stan")
source_file = os.path.join(path_to_this_file, "data", "sourcedata.h5")
uhecr_file = os.path.join(path_to_this_file, "data", "UHECRdata.h5")
table_path = os.path.join(path_to_this_file, "tables")
output_path = os.path.join(path_to_this_file, "output")

# make output path if it doesnt exist
if not os.path.exists(output_path):
    os.mkdir(output_path)

parser = argparse.ArgumentParser(description='Fit Stan Model to Data')
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
    help='The stan model considered (from arrival_direction, joint, joint_gmf)',
    choices=["arrival_direction", "joint", "joint_gmf"])
parser.add_argument(
    '--dtype',
    dest='dtype',
    action='store',
    default="real",
    type=str,
    help="Fit with simulated or real data (choose between 'sim' and 'real')",
    choices=["sim", "real"])
parser.add_argument(
    '--ptype',
    dest='ptype',
    action='store',
    default="p",
    type=str,
    help=
    "Type of particle used for back propagation (only used with joint_gmf model)."
)
parser.add_argument("--seed",
                    dest="seed",
                    action="store",
                    default=19990308,
                    type=int,
                    help="Random seed used for MCMC in Stan.")

parser.add_argument(
    '--sim_model',
    dest='sim_modeltype',
    action='store',
    default=None,
    help='The simulation model considered (from joint, joint_gmf)',
    choices=["joint", "joint_gmf"])
parser.add_argument("--tight_B",
                    dest="tight_B",
                    action="store_true",
                    help="Enable model with tighetened EGMF.",
                    default=False)


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


if __name__ == "__main__":

    args = parser.parse_args()

    # get filenames
    table_file = os.path.join(
        table_path, 'tables_{0}_{1}.h5'.format(args.source_type,
                                               args.detector_type))

    if args.sim_modeltype is not None:

        if args.tight_B:
            analysis_output_file = os.path.join(
                output_path,
                "{0}_fit_{5}_{1}_{2}_{3}_{4}_{6}_tightB.h5".format(
                    args.model_type, args.source_type, args.detector_type,
                    args.seed, args.ptype, args.dtype, args.sim_modeltype))
        else:
            analysis_output_file = os.path.join(
                output_path, "{0}_fit_{5}_{1}_{2}_{3}_{4}_{6}.h5".format(
                    args.model_type, args.source_type, args.detector_type,
                    args.seed, args.ptype, args.dtype, args.sim_modeltype))
    else:
        if args.tight_B:
            analysis_output_file = os.path.join(
                output_path,
                "tmp_{0}_fit_{5}_{1}_{2}_{3}_{4}_tightB.h5".format(
                    args.model_type, args.source_type, args.detector_type,
                    args.seed, args.ptype, args.dtype))
        else:
            analysis_output_file = os.path.join(
                output_path, "tmp_{0}_fit_{5}_{1}_{2}_{3}_{4}.h5".format(
                    args.model_type, args.source_type, args.detector_type,
                    args.seed, args.ptype, args.dtype))

    # get things related to detector
    detector_properties, alpha_T, M, Eth = get_detectorimports(
        args.detector_type)

    # create Data object
    data = Data()

    # add things to Data, method depends on sim vs real data
    if args.dtype == "sim":  # get data from sim_output_file
        # simulated data uses joint model, unless considering gmf propagation too
        if args.sim_modeltype == "joint_gmf":
            sim_output_file = os.path.join(
                output_path,
                "{0}_sim_{1}_{2}_{3}_{4}.h5".format(args.sim_modeltype,
                                                    args.source_type,
                                                    args.detector_type,
                                                    args.seed, args.ptype))
        # else:
        #     sim_output_file = os.path.join(
        #         output_path,
        #         "{0}_sim_{1}_{2}_{3}.h5".format(args.sim_modeltype,
        #                                         args.source_type,
        #                                         args.detector_type, args.seed))
        data.from_file(sim_output_file)

    elif args.dtype == "real":  # add source / UHECR / detector data manually
        data.add_source(source_file, args.source_type)
        data.add_uhecr(uhecr_file, args.detector_type, args.ptype)
        data.add_detector(detector_properties)

    # create Model, compile it, and set input
    if args.tight_B:
        model_name = os.path.join(
            stan_path,
            '{0}_model_composition_tightB.stan'.format(args.model_type))

    else:
        model_name = os.path.join(
            stan_path, '{0}_model_composition.stan'.format(args.model_type))

    print("tightB: ", args.tight_B, model_name)

    model = Model(model_filename=model_name, include_paths=stan_path)
    model.compile(reset=False)
    model.input(Eth=Eth)  # in EeV

    # perform the analysis
    summary = b'Fitting the model to given data.'
    analysis = Analysis(data,
                        model,
                        analysis_type=args.model_type,
                        filename=analysis_output_file,
                        summary=summary)

    # Each catalogue has a file of pre-computed values
    analysis.use_tables(table_file)

    # Fit the Stan model
    fit = analysis.fit_model(chains=16, iterations=500, seed=args.seed)

    # Save to analysis file
    analysis.save()
