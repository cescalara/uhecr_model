'''
Plot output figures from output files.
'''

import os
import argparse
import sys

from figures.output_figures import OutputFigures, SourceUHECRDist
from figures.deflection_figures import DeflectionFigures

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
    help=
    'either plot simulation results, fit sim results, fit data results, or do it all.',
    choices=['sim', 'fit_sim', 'fit_data', 'defl', 'all'])

parser.add_argument(
    '--dist',
    dest='dist_mode',
    action='store',
    default='all',
    help='Comparison type, compare with selected category in distribution.',
    choices=['model', 'source', 'detector', 'ptype', 'seed'])

parser.add_argument('--plot_all',
                    dest='plot_all',
                    action='store_true',
                    help='Run all, with all possible configurations.')

parser.add_argument('-d',
                    '--debug',
                    dest="debug",
                    action="store",
                    default=None,
                    help="Plot in debug mode with TA, SBG, proton, all models",
                    choices=["sim", "fit_sim", "fit_data", "all", None])

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

    if args.plot_all:  # use maximum configuration
        sources = ["SBG_23", "2FHL_250Mpc", "swift_BAT_213"]
        detectors = ["TA2015", "auger2014"]
        sim_models = ["joint", "joint_gmf"]
        sim_models_for_fit = ["joint", "joint_gmf"]
        fit_models = ["arrival_direction", "joint", "joint_gmf"]
        ptypes = ["p", "He", "N", "Si", "Fe"]
        seeds = [19990308, 4968460, 165490]
        end_labels = [None, "tight_B"]
        verbose = False
        header = "all"

    elif args.debug is not None:  # run the debug case
        sources = ["SBG_23"]
        detectors = ["TA2015"]
        sim_models = ["joint_gmf"]
        sim_models_for_fit = ["joint_gmf"]
        fit_models = ["arrival_direction", "joint", "joint_gmf"]
        ptypes = ["N"]
        seeds = [19990308]
        end_labels = [None]
        verbose = True
        header = "debug"
        args.mode = args.debug

    else:  # use config manually set up over here
        sources = ["SBG_23"]
        detectors = ["TA2015"]
        sim_models = ["joint"]
        sim_models_for_fit = ["joint"]
        fit_models = ["arrival_direction", "joint", "joint_gmf"]
        ptypes = ["p", "N"]
        # seeds = [19990308, 4968460, 165490]
        seeds = [19990308]
        end_labels = [None]
        verbose = args.verbose
        header = "noIGMF"

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
                                str(end_labels) + "\nHeader: " + header +
                                "\nDo you want to continue (y/(n))?\n")

            if confirm_run != "y":
                sys.exit("Terminating plot construction.")

    # setup configurations for GMF deflection plots
    gmf_model = "JF12"

    # setup the configuration list, which contains a list of tuples that is used for the simulation / fitting
    # order will be like: (source, detector, model, ptype, seed)
    if args.mode == "sim" or args.mode == "all":
        sim_config_list = [(src, d, m, p, se, end) for src in sources
                           for d in detectors for m in sim_models
                           for p in ptypes for se in seeds
                           for end in end_labels]
    if args.mode.find("fit") != -1 or args.mode == "all":
        fit_config_list = [(src, d, m, p, se, end) for src in sources
                           for d in detectors for m in fit_models
                           for p in ptypes for se in seeds
                           for end in end_labels]

    # set up the dictionary that contains the list of configurations.
    # this will be read to fdist class, and then one can create different distributions
    # with each comparison type
    # the list is setup from the config lists
    config_dict = {
        "model": fit_models,
        "source": sources,
        "detector": detectors,
        "ptype": ptypes,
        "seed": seeds,
        "end_label": end_labels,
    }

    # # run the actual thing
    # if args.mode == "sim" or args.mode == "all":
    #     for (source, detector, sim_model, ptype, seed,
    #          end_label) in sim_config_list:
    #         sim_args = {
    #             "source": source,
    #             "detector": detector,
    #             "model": sim_model_for_fit,
    #             "ptype": ptype,
    #             "seed": seed,
    #             "end_label": end_label
    #         }
    #         curr_config = "simulation: Source: {0}, Detector: {1}, sim_model: {2}, ptype: {3}, seed: {4}, additional config: {5}".format(
    #             source, detector, sim_model, ptype, seed, end_label)
    #         print("Current configuration:", curr_config)

    #         print("Plotting simulation results with: ", sim_model)
    #         output_fig_creator = OutputFigures(sim_args,
    #                                            sim_model=sim_model,
    #                                            verbose=verbose)
    #         if not args.dryrun:
    #             output_fig_creator.src_uhecr_skymap()
    #             output_fig_creator.corner()

    #         # TODO: if GMF sims, then have some way to plot the simulation figures

    if args.mode.find("fit") != -1 or args.mode == "all":
        for (source, detector, fit_model, ptype, seed,
             end_label) in fit_config_list:

            fig_args = {
                "source": source,
                "detector": detector,
                "model": fit_model,
                "ptype": ptype,
                "seed": seed,
                "end_label": end_label
            }

            if args.mode == "fit_sim" or args.mode == "fit" or args.mode == "all":
                for sim_model_for_fit in sim_models_for_fit:
                    print("Plotting results from simulated model: ",
                          sim_model_for_fit)
                    curr_config = "Fit simulation: Source: {0}, Detector: {1}, fit_model: {2}, ptype: {3}, seed: {4}, sim_model: {5}, additional config: {6}".format(
                        source, detector, fit_model, ptype, seed,
                        sim_model_for_fit, end_label)

                    print("Current configuration:", curr_config)

                    output_fig_creator = OutputFigures(
                        fig_args,
                        sim_model=sim_model_for_fit,
                        header=header,
                        verbose=verbose)

                    src_dist_creator = SourceUHECRDist(
                        fig_args,
                        config_dict,
                        sim_model=sim_model_for_fit,
                        header=header,
                        verbose=verbose)

                    if not args.dryrun:

                        output_fig_creator.src_uhecr_skymap()
                        output_fig_creator.corner()

                        if fit_model.find("joint") != -1:
                            output_fig_creator.dist(param="B")
                            output_fig_creator.association_skymap()

                        src_dist_creator.plotdist(compare_label="model",
                                                  extend=True)

                        src_dist_creator.plotdist(compare_label="model",
                                                  cumul=True,
                                                  extend=True)

            elif args.mode == "fit_data" or args.mode == "fit" or args.mode == "all":
                print("Plotting results with data from: ", detector)
                curr_config = "Fit data: Source: {0}, Detector: {1}, fit_model: {2}, ptype: {3}, seed: {4}, additional config: {5}".format(
                    source, detector, fit_model, ptype, seed, end_label)

                print("Current configuration:", curr_config)

                output_fig_creator = OutputFigures(fig_args,
                                                   header=header,
                                                   verbose=verbose)

                src_dist_creator = SourceUHECRDist(fig_args,
                                                   config_dict,
                                                   header=header,
                                                   verbose=verbose)
                if not args.dryrun:
                    output_fig_creator.src_uhecr_skymap()

                    if fit_model.find("joint") != -1:
                        output_fig_creator.association_skymap()
                        output_fig_creator.dist(param="f")

                    output_fig_creator.corner()

                    src_dist_creator.plotdist(compare_label="model",
                                              extend=True)

            print("\n")

    # create deflection figures
    if args.mode == "defl" or args.mode == "all":

        for detector in detectors:
            defl_fig_creator = DeflectionFigures(detector=detector,
                                                 header=header,
                                                 gmf_model=gmf_model)

            defl_fig_creator.plot_smearing(ptype="p",
                                           nodefl=True)  # no deflections
            for ptype in ptypes:
                defl_fig_creator.plot_smearing(ptype=ptype)
                defl_fig_creator.plot_smearing(ptype=ptype, sel_uhecr_idx=24)

                defl_fig_creator.plot_circles(ptype=ptype)
                defl_fig_creator.plot_circles(ptype=ptype, sel_uhecr_idx=24)