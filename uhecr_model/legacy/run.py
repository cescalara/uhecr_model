'''Run different scripts based on different configurations etc.'''

import os
import argparse
from tabnanny import verbose

from simulate_data_new import Simulation
from fit_model_new import FitModel
from figures.output_figures_new import OutputFigures, SourceUHECRDist

# paths to important files
path_to_this_file = os.path.abspath(os.path.dirname(__file__))
plot_path = os.path.join(path_to_this_file, "..", "..", "plots")
result_path = os.path.join(path_to_this_file, "..", "..", "results")

# make corresponding directories if they dont exist
if not os.path.exists(result_path):
    os.mkdir(result_path)

parser = argparse.ArgumentParser(
    description='Run (Simulation + ) Fit with the models.')

parser.add_argument(
    '--mode',
    dest='mode',
    action='store',
    default='fit_data',
    help='either simulate, fit with simulation, fit with data, or do it all.',
    choices=['sim', 'fit_sim', 'fit_data', 'fit', 'all'])

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
    choices=["sim", "fit_sim", "fit_data", "plot", "all", None])

parser.add_argument(
    '--dryrun',
    dest="dryrun",
    action="store_true",
    help="Perform dryrun and show what would be simulated / fitted.")

if __name__ == "__main__":
    '''
    list of possible configurations
    to modify, need to manually modify the list

    sources = ["SBG_23", "2FHL_250Mpc", "swift_BAT_213"]
    detectors = ["TA2015", "auger2014"]
    sim_models = ["joint", "joint_gmf"]
    fit_models = ["arrival_direction", "joint", "joint_gmf"]
    ptypes = ["p", "He", "N", "Si", "Fe"]
    seeds = [19990308, 4968460, 165490]
    '''

    args = parser.parse_args()

    # currently i can only think of setting this manually, maybe theres a better way?
    sources = ["2FHL_250Mpc", "SBG_23"]
    detectors = ["TA2015", "auger2014"]
    sim_models = ["joint_gmf"]
    fit_models = ["arrival_direction", "joint", "joint_gmf"]
    ptypes = ["p", "N"]
    seeds = [19990308, 4968460, 165490]

    # text file that shows the runs that have been completed. This is way better than looking through
    # stdout and stderr.
    progress_file = os.path.join(result_path, "progress_report.txt")

    # if it exists, remove it to reset it
    if os.path.exists(progress_file):
        os.remove(progress_file)

    # manually setting the simulation used to fit. There is yet again a better way to do this.
    sim_model_for_fit = "joint_gmf"
    # sim_model_for_fit = sim_models[0]

    # some label at the end if we want to perform additional modifications to the process.
    # set this to whatever suits your needs.
    # make sure to modify whatever u need to modify within the code as well!
    end_labels = [None, "tight_B"]

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

    # dryrun, simply print the output and dont run anything
    if args.dryrun:
        if args.mode == "sim" or args.mode == "all":
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

                curr_config = "simulation: Source: {0}, Detector: {1}, sim_model: {2}, ptype: {3}, seed: {4}".format(
                    source, detector, sim_model, ptype, seed)
                print("Current configuration:", curr_config)

                print("Simulating the model: ", sim_model)
                simulation = Simulation(sim_args, verbose=True)
                print("\n")

                with open(progress_file, "a") as f:
                    f.write("\n")
                    f.write(curr_config)

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

                if args.mode == "fit_sim" or args.mode == "fit" or args.mode == "all":
                    print("Fitting with simulated model: ", sim_model_for_fit)
                    curr_config = "Fit simulation: Source: {0}, Detector: {1}, fit_model: {2}, ptype: {3}, seed: {4}, sim_model: {5}".format(
                        source, detector, fit_model, ptype, seed,
                        sim_model_for_fit)

                    print("Current configuration:", curr_config)

                    fit = FitModel(fit_args,
                                   sim_model=sim_model_for_fit,
                                   verbose=True)

                    with open(progress_file, "a") as f:
                        f.write("\n")
                        f.write(curr_config)

                print("\n")

                if args.mode == "fit_data" or args.mode == "fit" or args.mode == "all":
                    print("Fitting with data")
                    curr_config = "Fit data: Source: {0}, Detector: {1}, fit_model: {2}, ptype: {3}, seed: {4}".format(
                        source, detector, fit_model, ptype, seed)

                    print("Current configuration:", curr_config)

                    fit = FitModel(fit_args, verbose=True)

                    with open(progress_file, "a") as f:
                        f.write("\n\n")
                        f.write(curr_config)

                print("\n")

    # run a debug run with TA, SBG, protons, one seed only (for data)
    # for simulation, run with TA, SBG, protons, multiple seeds, joint model
    elif args.debug is not None:
        seed_list = [19990308, 4968460, 165490]
        for seed in seeds:

            if args.debug == "sim" or args.debug == "all":
                sim_args = {
                    "source": "SBG_23",
                    "detector": "TA2015",
                    "model": "joint",
                    "ptype": "p",
                    "seed": seed,
                    "end_label": None
                }
                print("Current configuration:")
                print(
                    "Source: {0}, Detector: {1}, sim_model: {2}, ptype: {3}, seed: {4}"
                    .format("SBG_23", "TA2015", "joint", "p", seed))

                print("Simulating the model: ", "joint")
                simulation = Simulation(sim_args, verbose=True)
                simulation.set_simulation()
                simulation.simulate()
                simulation.save()

            if args.debug == "fit_sim" or args.debug == "all":
                fit_args = {
                    "source": "SBG_23",
                    "detector": "TA2015",
                    "model": "joint",
                    "ptype": "p",
                    "seed": seed,
                    "end_label": None
                }

                print("Fitting with simulated model: joint")
                fit_sim = FitModel(fit_args, sim_model="joint", verbose=True)
                fit_sim.set_analysis()
                fit_sim.fit()
                fit_sim.save()

        if args.debug == "fit_data" or args.debug == "all":
            fit_args = {
                "source": "SBG_23",
                "detector": "TA2015",
                "model": "joint",
                "ptype": "p",
                "seed": seed_list[0],
                "end_label": None
            }

            print("Fitting with data")
            fit_data = FitModel(fit_args, verbose=True)
            fit_data.set_analysis()
            fit_data.fit()
            fit_data.save()

        if args.debug == "plot" or args.debug == "all":
            # plot results
            print("Plotting results")

            fig_args = {
                "source": "SBG_23",
                "detector": "TA2015",
                "model": "joint",
                "ptype": "p",
                "seed": seed_list,
                "end_label": None
            }

            fit_args = {
                "source": "SBG_23",
                "detector": "TA2015",
                "model": "joint",
                "ptype": "p",
                "seed": seed_list[0],
                "end_label": None
            }
            # simulation
            print("Plotting simulation...")
            output_fig_creator = OutputFigures(fit_args,
                                               sim_model="joint",
                                               verbose=True)

            output_fig_creator.src_uhecr_skymap()
            output_fig_creator.corner()

            src_dist_creator = SourceUHECRDist(fig_args,
                                               sim_model="joint",
                                               verbose=True)

            src_dist_creator.plotdist(seed_list, extend=True)

            # data
            print("Plotting data...")
            output_fig_creator = OutputFigures(fit_args, verbose=True)

            output_fig_creator.src_uhecr_skymap()
            output_fig_creator.association_skymap()
            output_fig_creator.corner()

            src_dist_creator = SourceUHECRDist(fig_args, verbose=True)

            src_dist_creator.plotdist(seed_list, extend=True)

    # run the actual thing
    else:
        if args.mode == "sim" or args.mode == "all":
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
                curr_config = "simulation: Source: {0}, Detector: {1}, sim_model: {2}, ptype: {3}, seed: {4}".format(
                    source, detector, sim_model, ptype, seed)
                print("Current configuration:", curr_config)

                print("Simulating the model: ", sim_model)
                simulation = Simulation(sim_args, verbose=True)

                # modify B-field to those that converge with N
                B = 12 if end_label == "tight_B" else 20

                simulation.set_simulation(B=B)
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

                if fit_model.find("joint") != -1:
                    model_fname = "{0}_model_tightB.stan".format(
                        fit_model) if end_label == "tight_B" else None
                else:
                    model_fname = None

                print("Running: ", args.mode)
                print("Current configuration:")
                print(
                    "Source: {0}, Detector: {1}, fit_model: {2}, ptype: {3}, seed: {4}"
                    .format(source, detector, fit_model, ptype, seed))

                if args.mode == "fit_sim" or args.mode == "fit" or args.mode == "all":
                    print("Fitting with simulated model: ", sim_model_for_fit)
                    curr_config = "Fit simulation: Source: {0}, Detector: {1}, fit_model: {2}, ptype: {3}, seed: {4}, sim_model: {5}".format(
                        source, detector, fit_model, ptype, seed,
                        sim_model_for_fit)

                    print("Current configuration:", curr_config)

                    fit = FitModel(fit_args,
                                   sim_model=sim_model_for_fit,
                                   verbose=False)

                    fit.set_analysis(model_fname=model_fname)
                    fit.fit()
                    fit.save()

                    with open(progress_file, "a") as f:
                        f.write("\n")
                        f.write(curr_config)

                if args.mode == "fit_data" or args.mode == "fit" or args.mode == "all":
                    print("Fitting with data")
                    curr_config = "Fit data: Source: {0}, Detector: {1}, fit_model: {2}, ptype: {3}, seed: {4}".format(
                        source, detector, fit_model, ptype, seed)

                    print("Current configuration:", curr_config)

                    fit = FitModel(fit_args, verbose=False)

                    fit.set_analysis(model_fname=model_fname)
                    fit.fit()
                    fit.save()

                    with open(progress_file, "a") as f:
                        f.write("\n\n")
                        f.write(curr_config)

                print("\n")
