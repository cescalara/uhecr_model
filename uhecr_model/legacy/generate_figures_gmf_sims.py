'''
Python script to generate figures from output files.
Converted from figures.ipynb.
'''
import os
import argparse
import csv

from figures.output_figures import OutputFigures, SourceUHECRDist

path_to_this_file = os.path.abspath(os.path.dirname(__file__))
result_path = os.path.join(path_to_this_file, "..", "..", "result")
plot_path = os.path.join(path_to_this_file, "..", "..", "plots")
output_path = os.path.join(path_to_this_file, "output")

# create figure directory if it doesnt exist
if not os.path.exists(plot_path):
    os.mkdir(plot_path)

if not os.path.exists(result_path):
    os.mkdir(result_path)

if __name__ == "__main__":
    # for joint + gmf model, data comparisons only
    sources = ["SBG_23", "2FHL_250Mpc", "swift_BAT_213"]
    # sources = ["SBG_23"]
    # detectors = ["TA2015", "auger2014"]
    models = ["arrival_direction", "joint", "joint_gmf"]
    # dtypes = ["real", "sim"]
    # ptypes = ["p", "He", "N", "O", "Si", "Fe"]
    ptypes = ["p", "He", "N", "O"]
    # ptypes = ["p"]
    seeds = [
        19990308, 16852056, 65492186, 9999999, 9953497, 163591, 78520492,
        51610230, 8946064, 9526433, 84910385, 11223344, 2165490489, 7498365
    ]
    sim_model = "joint_gmf"

    output_file = "{0}_fit_{1}_{2}_{3}_{4}_{5}_{6}.h5"
    sim_output_file = "{0}_sim_{1}_{2}_{3}_{4}_{5}.h5"
    model_labels = ["Arrival Direction", "Joint", "Joint + GMF"]

    # # GMF sims: arrival vs joint vs joint_gmf, all sources, all particle types, all seeds, TA only
    # savefile_dist = "GMFplots_{0}_{1}_{2}_{3}_{4}_dist_{5}sim.png"
    # for ptype in ptypes:
    #     for source in sources:
    #         for seed in seeds:
    #             print("config: all models, source: {0}, seed: {1}, ptype: {2}".
    #                   format(source, seed, ptype))

    #             output_files = []
    #             for model in models:
    #                 if model == "joint_gmf":
    #                     output_files.append(
    #                         os.path.join(
    #                             output_path,
    #                             output_file.format(model, "sim", source,
    #                                                "TA2015", seed, ptype,
    #                                                sim_model)))
    #                 else:
    #                     output_files.append(
    #                         os.path.join(
    #                             output_path,
    #                             output_file.format(model, "sim", source,
    #                                                "TA2015", seed, "p",
    #                                                sim_model)))

    #             source_label = "Swift-BAT" if source == "swift_BAT_213" else source.split(
    #                 "_")[0]

    #             sim_output_file = os.path.join(
    #                 output_path,
    #                 sim_output_file.format(sim_model, source, "TA2015", seed,
    #                                        ptype, "notightB"))

    #             dist_creator = SourceUHECRDist(sim_output_file=sim_output_file)

    #             dist_creator.get_fs(output_files)

    #             dist_creator.plotdist(model_labels, title=source_label)
    #             dist_creator.save(
    #                 os.path.join(
    #                     plot_path,
    #                     savefile_dist.format("sim", source, "TA2015", seed,
    #                                          ptype, "notightB")))

    # # GMF sims: cumulative distribution
    # savefile_dist = "GMFplots_{0}_{1}_{2}_{3}_{4}_cumul_dist_{5}sim.png"
    # for ptype in ptypes:
    #     output_files = []
    #     for model in models:
    #         print("config: cumulative, model: {2}, source: {0}, ptype: {1}".
    #               format("SBG_23", ptype, model))
    #         output_files_per_model = []
    #         for seed in seeds:

    #             if model == "joint_gmf":
    #                 output_files_per_model.append(
    #                     os.path.join(
    #                         output_path,
    #                         output_file.format(model, "sim", "SBG_23",
    #                                            "TA2015", seed, ptype,
    #                                            sim_model)))
    #             else:
    #                 output_files_per_model.append(
    #                     os.path.join(
    #                         output_path,
    #                         output_file.format(model, "sim", "SBG_23",
    #                                            "TA2015", seed, "p",
    #                                            sim_model)))

    #         output_files.append(output_files_per_model)

    #     sim_output_file = os.path.join(
    #         output_path,
    #         sim_output_file.format(sim_model, "SBG_23", "TA2015", seed, ptype,
    #                                "notightB"))
    #     dist_creator = SourceUHECRDist(sim_output_file=sim_output_file)

    #     dist_creator.get_fs(output_files, cumul=True)

    #     dist_creator.plotdist(model_labels, title="SBG")
    #     dist_creator.save(
    #         os.path.join(
    #             plot_path,
    #             savefile_dist.format("sim", "SBG_23", "TA2015", ptype,
    #                                  sim_model, "notightB")))

    # # GMF sims: corner plots, all seeds, all sources, all particle types, TA only
    # savefile_corner = "GMFplots_{0}_{1}_{2}_{3}_{4}_{5}_{6}_corner_{7}sim.png"
    # for ptype in ptypes:
    #     for source in sources:
    #         for seed in seeds:
    #             for model in models:
    #                 print(
    #                     "corner plots: model: {3}, source: {0}, seed: {1}, ptype: {2}"
    #                     .format(source, seed, ptype, model))

    #                 sim_output_file = os.path.join(
    #                     output_path,
    #                     sim_output_file.format(sim_model, source, "TA2015",
    #                                            seed, ptype, "notightB"))

    #                 fig_args = {
    #                     "model": model,
    #                     "dtype": "sim",
    #                     "source": source,
    #                     "detector": "TA2015",
    #                     "seed": seed,
    #                     "ptype": ptype,
    #                 }

    #                 output_fig_creator = OutputFigures(
    #                     fig_args,
    #                     sim_model=sim_model,
    #                     sim_output_file=sim_output_file)

    #                 output_fig_creator.corner_sim(savefile=os.path.join(
    #                     plot_path,
    #                     savefile_corner.format(model, "sim", source, "TA2015",
    #                                            seed, ptype, sim_model,
    #                                            "notightB")))

    # sources = ["SBG_23", "2FHL_250Mpc", "swift_BAT_213"]
    sources = ["SBG_23"]
    # detectors = ["TA2015", "auger2014"]
    models = ["arrival_direction", "joint", "joint_gmf"]
    # dtypes = ["real", "sim"]
    # ptypes = ["p", "He", "N", "O", "Si", "Fe"]
    # ptypes = ["p", "He", "N", "O"]
    ptypes = ["p"]
    seeds = [
        19990308, 16852056, 65492186, 9999999, 9953497, 163591, 78520492,
        51610230, 8946064, 9526433, 84910385, 11223344, 7498365, 196560486,
        98765432, 318979
    ]
    # seeds = [19990308, 16852056, 65492186, 9999999, 9953497]
    sim_model = "joint_gmf"

    output_file = "{0}_fit_{1}_{2}_{3}_{4}_{5}_{6}.h5"
    sim_output_file = "{0}_sim_{1}_{2}_{3}_{4}_{5}.h5"
    model_labels = ["Arrival Direction", "Joint", "Joint + GMF"]

    # GMF sims: cumulative distribution
    savefile_dist = "GMFplots_{0}_{1}_{2}_{3}_{4}_cumul_dist_{5}sim.png"
    for ptype in ptypes:
        output_files = []
        for model in models:
            print("config: cumulative, model: {2}, source: {0}, ptype: {1}".
                  format("SBG_23", ptype, model))
            output_files_per_model = []
            for seed in seeds:

                if model == "joint_gmf":
                    output_files_per_model.append(
                        os.path.join(
                            output_path,
                            output_file.format(model, "sim", "SBG_23",
                                               "TA2015", seed, ptype,
                                               sim_model)))
                else:
                    output_files_per_model.append(
                        os.path.join(
                            output_path,
                            output_file.format(model, "sim", "SBG_23",
                                               "TA2015", seed, "p",
                                               sim_model)))

            output_files.append(output_files_per_model)

        sim_output_file = os.path.join(
            output_path,
            sim_output_file.format(sim_model, "SBG_23", "TA2015", seed, ptype,
                                   "notightB"))
        dist_creator = SourceUHECRDist(sim_output_file=sim_output_file)

        dist_creator.get_fs(output_files, cumul=True)

        dist_creator.plotdist(model_labels, title="SBG")
        dist_creator.save(
            os.path.join(
                plot_path,
                savefile_dist.format("sim", "SBG_23", "TA2015", ptype,
                                     sim_model, "notightB")))

    # GMF sims: corner plots, cumulative, all sources, all particle types, TA only
    savefile_corner_cumul = "GMFplots_{0}_{1}_{2}_{3}_{4}_{5}_corner_cumul_{6}sim.png"
    for ptype in ptypes:
        for source in sources:
            for model in ["joint_gmf"]:
                print(
                    "corner plots, cumul: model: {2}, source: {0}, ptype: {1}".
                    format(source, ptype, model))

                sim_output_file = os.path.join(
                    output_path,
                    sim_output_file.format(sim_model, source, "TA2015",
                                           19990308, ptype, "notightB"))

                fig_args = {
                    "model": model,
                    "dtype": "sim",
                    "source": source,
                    "detector": "TA2015",
                    "seed": 1990308,  # can be anything since this is ignored
                    "ptype": ptype
                }

                output_fig_creator = OutputFigures(
                    fig_args,
                    sim_model=sim_model,
                    sim_output_file=sim_output_file)

                output_file_list = []
                for seed in seeds:
                    output_file_list.append(
                        os.path.join(
                            output_path,
                            "{0}_fit_{5}_{1}_{2}_{3}_{4}_{6}.h5".format(
                                model, source, "TA2015", seed, ptype, "sim",
                                sim_model)))

                output_fig_creator.corner_sim(savefile=os.path.join(
                    plot_path,
                    savefile_corner_cumul.format(model, "sim", source,
                                                 "TA2015", ptype, sim_model,
                                                 "notightB")),
                                              cumul=True,
                                              output_files=output_file_list)
