'''
Python script to construct tables for each source and UHECR dataset.
Converted from precompute_tables.ipynb.
'''

import os
import h5py
import numpy as np
import argparse

from fancy import Data, Model, Analysis

parser = argparse.ArgumentParser(description='Precompute exposure integral and energy grid.')
parser.add_argument('--source', dest='source_type', action='store', default="SBG_23", type=str,
                    help='The source catalogue used (from SBG_23, 2FHL_250Mpc, swift_BAT_213)',
                    choices=["SBG_23", "2FHL_250Mpc", "swift_BAT_213"])
parser.add_argument('--detector', dest='detector_type', action='store', default="TA2015", type=str,
                    help='The type of detector config (from TA2015, auger2014, auger2010)',
                    choices=["TA2015", "auger2014", "auger2010"])
parser.add_argument('--model', dest='model_type', action='store', default="joint", type=str,
                    help='The stan model considered for simulation (from joint, joint_gmf)',
                    choices=["joint", "joint_gmf"])
parser.add_argument("--check", dest="check", action="store_true", 
                    help="Check if table is correctly constructed.")

def get_detectorimports(detector_type):
    '''Get variables imported by (detector_name).py'''
    if detector_type == "TA2015":
        from fancy.detector.TA2015 import detector_properties, Eth
    elif detector_type == "auger2014":
        from fancy.detector.auger2014 import detector_properties, Eth
    elif detector_type == "auger2010":
        from fancy.detector.auger2010 import detector_properties, Eth
    else:
        raise Exception("Undefined detector type!")

    return detector_properties, Eth

# paths to important files
path_to_this_file = os.path.abspath(os.path.dirname(__file__))
stan_path = os.path.join(path_to_this_file, "..", "stan")
source_file = os.path.join(path_to_this_file, "..", "data", "sourcedata.h5")
uhecr_file = os.path.join(path_to_this_file, "..", "data", "UHECRdata.h5")
# table_path = os.path.join(path_to_this_file, "tables")

if __name__ == "__main__":

    args = parser.parse_args()

    # filename for table
    table_file = os.path.join(path_to_this_file, "tables_{0}_{1}.h5".format(args.source_type, args.detector_type))

    # get things related to detector
    detector_properties, Eth = get_detectorimports(
        args.detector_type)


    # create Data object
    data = Data()
    data.add_source(source_file, args.source_type) 
    data.add_uhecr(uhecr_file, args.detector_type) 
    data.add_detector(detector_properties)

    # create Model object and set input
    model_name = os.path.join(
        stan_path, '{0}_model.stan'.format(args.model_type))  # not used so can be anything
    model = Model(model_filename = model_name, include_paths = stan_path)
    model.input(Eth = Eth) # EeV

    # create ANalysis object and build tables
    summary = b'Precomputation for tables'
    analysis = Analysis(data, model, analysis_type = 'joint', 
                        filename = table_file, summary = summary)

    # exposure integral
    analysis.build_tables(fit_only = True)
    analysis.tables.save(table_file)

    # grid for arrival energies
    analysis.build_energy_table(table_file = table_file)

    # check keys as verification if enabled
    if args.check:
        with h5py.File(table_file, "r") as f:
            print("Keys for table_file: ", f.keys())

            print("Exposure integral keys:")
            print(f["main"].keys())   # main exposure integral table

            print("Energy Grid keys:")
            print(f["energy"].keys())
