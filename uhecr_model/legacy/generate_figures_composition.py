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
    # sources = ["SBG_23", "2FHL_250Mpc", "swift_BAT_213"]
    sources = ["2FHL_250Mpc"]
    # sources = ["swift_BAT_213"]
    # detector = "auger2014"
    detector = "TA2015"
    # detectors = ["TA2015", "auger2014"]
    models = ["arrival_direction", "joint", "joint_gmf"]
    # models = ["arrival_direction", "joint"]
    # dtypes = ["real", "sim"]
    # ptypes = ["p", "He", "N", "O", "Si", "Fe"]
    # ptypes = ["Fe"]
    ptypes = ["p", "N"]
    # seeds = [19990308, 16852056]
    seeds = [19990308]

    output_file = "tmp_{0}_fit_{1}_{2}_{3}_{4}_{5}_tightB.h5"
    model_labels = ["Arrival Direction", "Joint", "Joint + GMF"]
    # model_labels = ["Arrival Direction", "Joint"]
    # GMF data: joint vs joint_gmf, all sources, all particle types, only 2 seeds, TA only
    savefile_dist = "tmp_Compoplots_{0}_{1}_{2}_{3}_{4}_dist.png"
    for ptype in ptypes:
        for source in sources:
            for seed in seeds:
                print(
                    "config: TA, all models, source: {0}, seed: {1}, ptype: {2}"
                    .format(source, seed, ptype))

                output_files = []
                for model in models:
                    if model == "joint_gmf":
                        output_files.append(
                            os.path.join(
                                output_path,
                                output_file.format(model, "real", source,
                                                   detector, seed, ptype)))
                    elif model == "joint":
                        output_files.append(
                            os.path.join(
                                output_path,
                                output_file.format(model, "real", source,
                                                   detector, seed, ptype)))
                    else:
                        output_files.append(
                            os.path.join(
                                output_path,
                                "{0}_fit_{1}_{2}_{3}_{4}_{5}.h5".format(
                                    model, "real", source, detector, seed,
                                    "p")))

                source_label = "Swift-BAT" if source == "swift_BAT_213" else source.split(
                    "_")[0]

                dist_creator = SourceUHECRDist()

                dist_creator.get_fs(output_files)

                dist_creator.plotdist(model_labels, title=source_label)
                dist_creator.save(
                    os.path.join(
                        plot_path,
                        savefile_dist.format("real", source, detector, seed,
                                             ptype)))

    # plot assossciation skymap + corner plots, all sources, all particle types, only 2 seeds, TA only

    result_args_list = []
    for model in ["joint", "joint_gmf"]:
        for ptype in ptypes:
            for source in sources:
                for seed in seeds:

                    fig_args = {
                        "model": model,
                        "dtype": "real",
                        "source": source,
                        "detector": detector,
                        "seed": seed,
                        "ptype": ptype,
                    }
                    print("config: " + "    ".join([
                        "{0} : {1}".format(key, value)
                        for key, value in fig_args.items()
                    ]))

                    output_figs = OutputFigures(fig_args, tight_B=True)

                    savefile_assos_skymap = os.path.join(
                        plot_path,
                        "Compoplots_{0}_{1}_{2}_{3}_{4}_{5}_assos_skymap.png".
                        format(*fig_args.values()))

                    savefile_corner = os.path.join(
                        plot_path,
                        "Compoplots_{0}_{1}_{2}_{3}_{4}_{5}_corner.png".format(
                            *fig_args.values()))

                    dominant_sources_dict = output_figs.association_skymap(
                        savefile_assos_skymap)

                    output_figs.corner_data(savefile_corner)

                    result_args = {
                        "model": model,
                        "dtype": "real",
                        "source": source,
                        "detector": detector,
                        "seed": seed,
                        "ptype": ptype,
                        "tight_B": True,
                        "dominant_sources": list(dominant_sources_dict.keys()),
                        "counts": list(dominant_sources_dict.values())
                    }

                    result_args_list.append(result_args)

    # # write dominant sources to some csv file
    # field_names = [
    #     "model", "dtype", "source", "detector", "seed", "ptype", "tight_B",
    #     "dominant_sources", "counts"
    # ]

    # with open(os.path.join(result_path, "dominant_source_list.csv"), "a") as f:
    #     writer = csv.DictWriter(f, fieldnames=field_names)
    #     for result_args in result_args_list:
    #         writer.writerow(result_args)

    # # for joint + gmf model, data comparisons only
    # sources = ["SBG_23", "2FHL_250Mpc", "swift_BAT_213"]
    # # detectors = ["TA2015", "auger2014"]
    # models = ["joint_gmf"]
    # # dtypes = ["real", "sim"]
    # ptypes = ["p", "He", "N", "O", "Si", "Fe"]
    # seeds = [19990308]

    # # Composition: source comparison, all detectors, all seeds
    # savefile_dist = "Compoplots_{0}_{1}_{2}_{3}_{4}_allsrc_dist.png"
    # src_labels = ["SBG", "2FHL", "Swift-BAT"]
    # detector_labels = ["TA"]
    # for i, model in enumerate(models):
    #     for ptype in ptypes:
    #         for seed in seeds:
    #             print("config: all sources, ptype : {0}, seed: {1}".format(
    #                 ptype, seed))

    #             output_files = []
    #             for source in sources:
    #                 if model == "joint_gmf":
    #                     output_files.append(
    #                         os.path.join(
    #                             output_path,
    #                             output_file.format(model, "real", source,
    #                                                "TA2015", seed, ptype)))
    #                 elif model == "joint":
    #                     output_files.append(
    #                         os.path.join(
    #                             output_path,
    #                             output_file.format(model, "real", source,
    #                                                "TA2015", seed, "p")))
    #                 else:
    #                     output_files.append(
    #                         os.path.join(
    #                             output_path,
    #                             "{0}_fit_{1}_{2}_{3}_{4}_{5}.h5".format(
    #                                 model, "real", source, "TA2015", seed,
    #                                 "p")))

    #             dist_creator = SourceUHECRDist()

    #             dist_creator.get_fs(output_files)

    #             dist_creator.plotdist(src_labels,
    #                                   title=model_labels[-1] + ", " + ptype,
    #                                   extend=True)
    #             dist_creator.save(
    #                 os.path.join(
    #                     plot_path,
    #                     savefile_dist.format(model, "real", "TA2015", seed,
    #                                          ptype)))

    # # GMF sims: joint vs joint_gmf, all sources, all particle types, only 2 seeds, TA only
    # sim_model = "joint_gmf"
    # output_file = "{0}_fit_{1}_{2}_{3}_{4}_{5}_tightB.h5"
    # savefile_dist = "Compoplots_{0}_{1}_{2}_{3}_{4}_dist.png"

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
    #                 sim_output_file.format(sim_model, "SBG_23", "TA2015", seed,
    #                                        ptype, "tightB"))

    #             dist_creator = SourceUHECRDist(sim_output_file=sim_output_file)

    #             dist_creator.get_fs(output_files)

    #             dist_creator.plotdist(model_labels, title=source_label)
    #             dist_creator.save(
    #                 os.path.join(
    #                     plot_path,
    #                     savefile_dist.format("sim", source, "TA2015", seed,
    #                                          ptype)))

    # # GMF sims: cumulative distribution
    # savefile_dist = "Compoplots_{0}_{1}_{2}_{3}_{4}_cumul_dist.png"
    # for ptype in ptypes:
    #     output_files = []
    #     for model in models:
    #         output_files_per_model = []
    #         for seed in seeds:
    #             print("config: cumulative, source: {0}, seed: {1}, ptype: {2}".
    #                   format("SBG_23", seed, ptype))

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
    #                                "tightB"))
    #     dist_creator = SourceUHECRDist(sim_output_file=sim_output_file)

    #     dist_creator.get_fs(output_files, cumul=True)

    #     dist_creator.plotdist(model_labels, title="SBG")
    #     dist_creator.save(
    #         os.path.join(
    #             plot_path,
    #             savefile_dist.format("sim", "SBG_23", "TA2015", ptype,
    #                                  sim_model)))
# '''
# Python script to generate figures from output files.
# Converted from figures.ipynb.
# '''
# import os
# import argparse
# import csv

# from figures.output_figures import OutputFigures, SourceUHECRDist

# path_to_this_file = os.path.abspath(os.path.dirname(__file__))
# result_path = os.path.join(path_to_this_file, "..", "..", "result")
# plot_path = os.path.join(path_to_this_file, "..", "..", "plots")
# output_path = os.path.join(path_to_this_file, "output")

# # create figure directory if it doesnt exist
# if not os.path.exists(plot_path):
#     os.mkdir(plot_path)

# if not os.path.exists(result_path):
#     os.mkdir(result_path)

# if __name__ == "__main__":
#     # for joint + gmf model, data comparisons only
#     sources = ["SBG_23", "2FHL_250Mpc", "swift_BAT_213"]
#     # detectors = ["TA2015", "auger2014"]
#     models = ["arrival_direction", "joint", "joint_gmf"]
#     # dtypes = ["real", "sim"]
#     ptypes = ["p", "He", "N", "O", "Si", "Fe"]
#     seeds = [19990308, 16852056]

#     output_file = "{0}_fit_{1}_{2}_{3}_{4}_{5}_tightB.h5"
#     model_labels = ["Arrival Direction", "Joint", "Joint + GMF"]

#     # GMF data: joint vs joint_gmf, all sources, all particle types, only 2 seeds, TA only
#     savefile_dist = "Compoplots_{0}_{1}_{2}_{3}_{4}_dist.png"

#     output_files = []
#     for model in models:
#         if model == "joint_gmf":
#             output_files.append(
#                 os.path.join(
#                     output_path,
#                     output_file.format(model, "real", "SBG_23", "TA2015",
#                                        19990308, "p")))
#         elif model == "joint":
#             output_files.append(
#                 os.path.join(
#                     output_path,
#                     output_file.format(model, "real", "SBG_23", "TA2015",
#                                        19990308, "p")))
#         else:
#             output_files.append(
#                 os.path.join(
#                     output_path, "{0}_fit_{1}_{2}_{3}_{4}_{5}.h5".format(
#                         model, "real", "SBG_23", "TA2015", 19990308, "p")))

#     dist_creator = SourceUHECRDist()

#     dist_creator.get_fs(output_files)

#     dist_creator.plotdist(model_labels, title="SBG")
#     dist_creator.save(
#         os.path.join(
#             plot_path,
#             savefile_dist.format("real", "SBG_23", "TA2015", 19990308, "p")))

#     # plot assossciation skymap + corner plots, all sources, all particle types, only 2 seeds, TA only

#     fig_args = {
#         "model": "joint_gmf",
#         "dtype": "real",
#         "source": "SBG_23",
#         "detector": "TA2015",
#         "seed": 19990308,
#         "ptype": "p"
#     }
#     print("config: " + "    ".join(
#         ["{0} : {1}".format(key, value) for key, value in fig_args.items()]))

#     output_figs = OutputFigures(fig_args, tight_B=True)

#     savefile_assos_skymap = os.path.join(
#         plot_path,
#         "Compoplots_{0}_{1}_{2}_{3}_{4}_{5}_assos_skymap.png".format(
#             *fig_args.values()))

#     savefile_corner = os.path.join(
#         plot_path, "Compoplots_{0}_{1}_{2}_{3}_{4}_{5}_corner.png".format(
#             *fig_args.values()))

#     dominant_sources_dict = output_figs.association_skymap(
#         savefile_assos_skymap)

#     output_figs.corner_data(savefile_corner)

#     result_args = {
#         "model": "joint_gmf",
#         "dtype": "real",
#         "source": "SBG_23",
#         "detector": "TA2015",
#         "seed": 19990308,
#         "ptype": "p",
#         "tight_B": True,
#         "dominant_sources": list(dominant_sources_dict.keys()),
#         "counts": list(dominant_sources_dict.values())
#     }
# # write dominant sources to some csv file
# field_names = [
#     "model", "dtype", "source", "detector", "seed", "ptype",
#     "dominant_sources", "counts"
# ]

# with open(os.path.join(result_path, "dominant_source_list.csv"), "a") as f:
#     writer = csv.DictWriter(f, fieldnames=field_names)
#     for result_args in result_args_list:
#         writer.writerow(result_args)
