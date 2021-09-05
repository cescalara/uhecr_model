'''Run different scripts based on different configurations etc.'''

import os
import argparse
import sys

from simulate_data import Simulation
from fit_model import FitModel

# paths to important files
path_to_this_file = os.path.abspath(os.path.dirname(__file__))
plot_path = os.path.join(path_to_this_file, "..", "..", "plots")
outputs_path = os.path.join(path_to_this_file, "outputs")

parser = argparse.ArgumentParser(
    description='Run (Simulation + ) Fit with the models.')

parser.add_argument(
    '--mode',
    dest='mode',
    action='store',
    default='all',
    help='either simulate, fit with simulation, fit with data, or do it all.',
    choices=['sim', 'fit_sim', 'fit_data', 'fit', 'sim_and_fit', 'all'])

parser.add_argument('--run_all',
                    dest='run_all',
                    action='store_true',
                    help='Run all, with all possible configurations.')

parser.add_argument(
    '-d',
    '--debug',
    dest="debug",
    action="store",
    default=None,
    help="Perform debug run with TA, SBG, proton, all models",
    choices=["sim", "sim_and_fit", "fit_sim", "fit_data", "fit", "all", None])

parser.add_argument(
    '-y',
    '--yes',
    dest="skip_confirm",
    action="store_true",
    help="Skip confirmation to check if the right configuration is performed.")

parser.add_argument(
    '-v',
    '--verbose',
    dest="verbose",
    action="store_true",
    help="Show lots of print out statements, usually for debug purposes.")

parser.add_argument(
    '--dryrun',
    dest="dryrun",
    action="store_true",
    help="Perform dryrun and show what would be simulated / fitted.")

if __name__ == "__main__":
    '''
    list of possible configurations
    to modify, need to manually modify the list

    sources = ["SBG_23", "2FHL_250Mpc", "swift_BAT_213"]  # source catalogue
    detectors = ["TA2015", "auger2014"]  # detector / UHECR dataset used
    sim_models = ["joint", "joint_gmf"]  # simulations to construct
    sim_models_for_fit = ["joint", "joint_gmf"]  # simulations to use for fitting
    fit_models = ["arrival_direction", "joint", "joint_gmf"]  # models to fit with
    ptypes = ["p", "He", "N", "Si", "Fe"]  # type of nuclei used
    seeds = [19990308, 4968460, 165490]  # seed
    '''

    args = parser.parse_args()

    # set up the configuration lists

    if args.run_all:  # use maximum configuration
        sources = ["SBG_23", "2FHL_250Mpc", "swift_BAT_213"]
        detectors = ["TA2015", "auger2014"]
        sim_models = ["joint", "joint_gmf"]
        sim_models_for_fit = ["joint", "joint_gmf"]
        fit_models = ["arrival_direction", "joint", "joint_gmf"]
        ptypes = ["p", "He", "N", "Si", "Fe"]
        seeds = [19990308, 4968460, 165490]
        end_labels = [None, "tight_B"]
        verbose = False

    elif args.debug is not None:  # run the debug case
        sources = ["2FHL_250Mpc"]
        detectors = ["auger2014"]
        sim_models = ["joint"]
        sim_models_for_fit = ["joint"]
        fit_models = ["arrival_direction", "joint", "joint_gmf"]
        ptypes = ["N"]
        seeds = [19990308]
        end_labels = ["tight_B"]
        verbose = True
        args.mode = args.debug

        # Nsim = 1000

    else:  # use config manually set up over here
        sources = ["SBG_23"]
        detectors = ["TA2015"]
        sim_models = ["joint"]
        sim_models_for_fit = ["joint"]
        fit_models = ["arrival_direction", "joint", "joint_gmf"]
        ptypes = ["p", "N"]
        # seeds = [
        #     19990308, 747072, 402309, 476859, 638121, 821056, 626445, 125326,
        #     568333, 606135, 214978, 858061, 41247, 556608, 508103, 716897,
        #     878370, 942554, 964306, 722605, 919183, 992879, 154451, 714282,
        #     437735, 519750, 390711
        # ]
        seeds = [19990308]
        end_labels = ["noIGMF"]
        verbose = args.verbose

        Nsim = 1250

        # when performing dryrun, enable verbosity
        if args.dryrun:
            verbose = True

        if args.skip_confirm:
            pass
        else:
            # show some confirmation message
            confirm_run = input("This is the current configuration:\n Mode: " +
                                args.mode + "\nSources:" + str(sources) +
                                "\n Detectors: " + str(detectors) +
                                "\n Models to Simulate: " + str(sim_models) +
                                "\n Simulations used for fitting: " +
                                str(sim_models_for_fit) +
                                "\n Models used for fitting: " +
                                str(fit_models) + "\n Particle types: " +
                                str(ptypes) + "\n seeds: " + str(seeds) +
                                "\n Additional configurations: " +
                                str(end_labels) +
                                "\nDo you want to continue (y/(n))?\n")

            if confirm_run != "y":
                sys.exit("Terminating run.")

    # text file that shows the runs that have been completed. This is way better than looking through
    # stdout and stderr.
    progress_file = os.path.join(outputs_path, "progress_report.txt")

    if args.debug:
        progress_file = os.path.join(outputs_path, "progress_report_debug.txt")

    # if it exists, remove it to reset it
    if os.path.exists(progress_file):
        os.remove(progress_file)

    # if os.path.exists(progress_file):
    #     progress_file = os.path.join(outputs_path, "progress_report_1.txt")

    # setup the configuration list, which contains a list of tuples that is used for the simulation / fitting
    # order will be like: (source, detector, model, ptype, seed)
    if args.mode == "sim" or args.mode == "sim_and_fit" or args.mode == "all":
        sim_config_list = [(src, d, m, p, se, end) for src in sources
                           for d in detectors for m in sim_models
                           for p in ptypes for se in seeds
                           for end in end_labels]
    if args.mode.find("fit") != -1 or args.mode == "all":
        fit_config_list = [(src, d, m, p, se, end) for src in sources
                           for d in detectors for m in fit_models
                           for p in ptypes for se in seeds
                           for end in end_labels]

    # run the actual thing
    if args.mode == "sim" or args.mode == "sim_and_fit" or args.mode == "all":
        for (source, detector, sim_model, ptype, seed,
             end_label) in sim_config_list:
            sim_args = {
                "source": source,
                "detector": detector,
                "model": sim_model,
                "ptype": ptype,
                "seed": seed,
                "end_label": end_label
            }
            curr_config = "simulation: Source: {0}, Detector: {1}, sim_model: {2}, ptype: {3}, seed: {4}, additional config: {5}".format(
                source, detector, sim_model, ptype, seed, end_label)
            print("Current configuration:", curr_config)

            print("Simulating the model: ", sim_model)
            simulation = Simulation(sim_args, verbose=verbose)

            # modify B-field to those that converge with N
            B = 3 if end_label == "tight_B" else 20

            simulation.set_simulation(B=B, Nsim=1000)

            if not args.dryrun:
                simulation.simulate()
                simulation.save()

            print("\n")

            with open(progress_file, "a") as f:
                f.write("\n")
                f.write(curr_config)

            print("\n")

    if args.mode.find("fit") != -1 or args.mode == "all":
        for (source, detector, fit_model, ptype, seed,
             end_label) in fit_config_list:

            fit_args = {
                "source": source,
                "detector": detector,
                "model": fit_model,
                "ptype": ptype,
                "seed": seed,
                "end_label": end_label
            }

            # setting model filename for tight Bfield priors
            model_fname = None
            if fit_model.find("joint") != -1:
                if end_label == "tight_B":
                    model_fname = "{0}_model_tightB.stan".format(fit_model)
                elif end_label == "noIGMF":
                    model_fname = "{0}_model_noIGMF.stan".format(fit_model)
                elif end_label == "limitL":
                    model_fname = "{0}_model_limitL.stan".format(fit_model)
                else:
                    model_fname = None

            if args.mode == "fit_sim" or args.mode == "sim_and_fit" \
                 or args.mode == "fit" or args.mode == "all":
                for sim_model_for_fit in sim_models_for_fit:
                    print("Fitting with simulated model: ", sim_model_for_fit)
                    curr_config = "Fit simulation: Source: {0}, Detector: {1}, fit_model: {2}, ptype: {3}, seed: {4}, sim_model: {5}, additional config: {6}".format(
                        source, detector, fit_model, ptype, seed,
                        sim_model_for_fit, end_label)

                    print("Current configuration:", curr_config)

                    fit = FitModel(fit_args,
                                   sim_model=sim_model_for_fit,
                                   verbose=verbose)

                    fit.set_analysis(model_fname=model_fname)

                    if not args.dryrun:
                        fit.fit()
                        fit.save()

                    with open(progress_file, "a") as f:
                        f.write("\n")
                        f.write(curr_config)

            if args.mode == "fit_data" or args.mode == "fit" or args.mode == "all":
                print("Fitting with data")
                curr_config = "Fit data: Source: {0}, Detector: {1}, fit_model: {2}, ptype: {3}, seed: {4}, additional config: {5}".format(
                    source, detector, fit_model, ptype, seed, end_label)

                print("Current configuration:", curr_config)

                fit = FitModel(fit_args, verbose=verbose)

                fit.set_analysis(model_fname=model_fname)

                if not args.dryrun:
                    fit.fit()
                    fit.save()

                with open(progress_file, "a") as f:
                    f.write("\n")
                    f.write(curr_config)
                    f.write("\n")

            print("\n")
