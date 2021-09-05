'''
Python script of the fitting process in the notebook run_simulation.ipynb.

This is made so that this can be run on command line as a bash script.
'''
import os
from chardet import detect
from fancy import Data, Model, Analysis

import argparse

# paths to important files
path_to_this_file = os.path.abspath(os.path.dirname(__file__))
stan_path = os.path.join(path_to_this_file, "stan")
source_file = os.path.join(path_to_this_file, "data", "sourcedata.h5")
uhecr_file = os.path.join(path_to_this_file, "data", "UHECRdata.h5")
table_path = os.path.join(path_to_this_file, "tables")
output_path = os.path.join(path_to_this_file, "outputs")
output_sim_path = os.path.join(output_path, "sim")
output_data_path = os.path.join(output_path, "data")

# make output path if it doesnt exist
if not os.path.exists(output_path):
    os.mkdir(output_path)

if not os.path.exists(output_data_path):
    os.mkdir(output_data_path)


class FitModel():
    def __init__(self, fit_args, sim_model=None, verbose=False):
        '''Fit model to either data or simulation.'''

        self.source = fit_args["source"]
        self.detector = fit_args["detector"]
        self.ptype = fit_args["ptype"]
        self.model = fit_args["model"]

        self.seed = fit_args["seed"]
        self.end_label = fit_args["end_label"]

        self.sim = False if sim_model is None else True
        self.sim_model = sim_model  # enable simulation or not

        self.verbose = verbose  # for print statements

        self._init_files()

    def _init_files(self):
        '''Initialize filenames for output files / table files'''
        self.table_file = os.path.join(
            table_path, 'tables_{0}_{1}.h5'.format(self.source, self.detector))

        if self.end_label is not None:
            analysis_output_file = "{0}_fit_{1}_{2}_{3}_{4}_{5}.h5".format(
                self.model, self.source, self.detector, self.seed, self.ptype,
                self.end_label)
        else:
            analysis_output_file = "{0}_fit_{1}_{2}_{3}_{4}.h5".format(
                self.model, self.source, self.detector, self.seed, self.ptype)

        if self.sim:
            if self.end_label is not None:
                sim_output_file = "{0}_{1}_{2}_{3}_{4}.h5".format(
                    self.source, self.detector, self.seed, self.ptype,
                    self.end_label)
            else:
                sim_output_file = "{0}_{1}_{2}_{3}.h5".format(
                    self.source, self.detector, self.seed, self.ptype)

            self.sim_output_file = os.path.join(output_sim_path,
                                                self.sim_model, "simulation",
                                                sim_output_file)

            # make required directories
            if not os.path.exists(
                    os.path.join(output_sim_path, self.sim_model, "fit")):
                os.mkdir(os.path.join(output_sim_path, self.sim_model, "fit"))

            self.analysis_output_file = os.path.join(output_sim_path,
                                                     self.sim_model, "fit",
                                                     analysis_output_file)

        else:
            self.analysis_output_file = os.path.join(output_data_path,
                                                     analysis_output_file)

        if self.verbose:
            print("Table output file: ", self.table_file)
            print("Analysis Output File: ", self.analysis_output_file)

            if self.sim:
                print("Reading from Simulation Output File: ",
                      self.sim_output_file)

    def _init_data_and_model(self, model_fname=None):
        '''Initialize Data object, which contains source / UHECR/ detector information'''

        detector_properties, Eth = self._get_detectorimports()

        if self.verbose:
            print("Energy threshold: ", Eth)

        data = Data()

        if self.sim:
            data.from_file(self.sim_output_file)

        else:
            data.add_source(source_file, self.source)
            data.add_uhecr(uhecr_file, self.detector)
            data.add_detector(detector_properties)

        model_fname = '{0}_model.stan'.format(
            self.model) if model_fname is None else model_fname

        model_filename = os.path.join(stan_path, model_fname)

        model = Model(model_filename=model_filename, include_paths=stan_path)
        model.compile(reset=False)
        model.input(Eth=Eth)  # in EeV

        if self.verbose:
            print("Model file used: ", model_filename)

        return data, model

    def _get_detectorimports(self):
        '''Get variables imported by (detector_name).py'''
        if self.detector == "TA2015":
            from fancy.detector.TA2015 import detector_properties, Eth
        elif self.detector == "auger2014":
            from fancy.detector.auger2014 import detector_properties, Eth
        elif self.detector == "auger2010":
            from fancy.detector.auger2010 import detector_properties, Eth
        else:
            raise Exception("Undefined detector type!")

        return detector_properties, Eth

    def set_analysis(self, model_fname=None):
        '''Setup the analysis which conducts the fitting'''

        data, model = self._init_data_and_model(
            model_fname=model_fname)  # initialize Data() object

        summary = b'Fitting the model to given data.'
        self.analysis = Analysis(data,
                                 model,
                                 analysis_type=self.model,
                                 filename=self.analysis_output_file,
                                 summary=summary)

        # Each catalogue has a file of pre-computed values
        self.analysis.use_tables(self.table_file)

    def fit(self, chains=24, iterations=400, warmup=100):
        '''Perform the fit.'''
        # Fit the Stan model
        fit = self.analysis.fit_model(chains=chains,
                                      iterations=iterations,
                                      warmup=warmup,
                                      seed=self.seed)

    def save(self):
        '''Save to output file'''
        # Save to analysis file
        self.analysis.save()
