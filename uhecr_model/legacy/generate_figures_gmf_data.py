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
    # detectors = ["TA2015", "auger2014"]
    detector = "TA2015"
    models = ["arrival_direction", "joint", "joint_gmf"]
    # dtypes = ["real", "sim"]
    # ptypes = ["p", "He", "N", "O", "Si", "Fe"]
    ptypes = ["p", "N"]
    # ptypes = ["p"]
    # seeds = [19990308, 16852056]
    seeds = [19990308]

    output_file = "tmp_{0}_fit_{1}_{2}_{3}_{4}_{5}.h5"
    model_labels = ["Arrival Direction", "Joint", "Joint + GMF"]
    # model_labels = ["Arrival Direction", "Joint"]

    # GMF data: joint vs joint_gmf, all sources, all particle types, only 2 seeds, TA only
    savefile_dist = "tmp_GMFplots_{0}_{1}_{2}_{3}_{4}_dist.png"
    for ptype in ptypes:
        for source in sources:
            for seed in seeds:
                print(
                    "config: TA, all models, source: {0}, seed: {1}, ptype: {2}"
                    .format(source, seed, ptype))

                output_files = []
                for model in models:

                    if model == "arrival_direction":
                        output_files.append(
                            os.path.join(
                                output_path,
                                "{0}_fit_{1}_{2}_{3}_{4}_{5}.h5".format(
                                    model, "real", source, detector, seed,
                                    "p")))

                    else:
                        output_files.append(
                            os.path.join(
                                output_path,
                                output_file.format(model, "real", source,
                                                   detector, seed, ptype)))

                # print(output_files)

                source_label = "Swift-BAT" if source == "swift_BAT_213" else source.split(
                    "_")[0]

                extend = True if source == "2FHL_250Mpc" else False

                dist_creator = SourceUHECRDist()

                dist_creator.get_fs(output_files)

                dist_creator.plotdist(model_labels,
                                      title=source_label,
                                      extend=extend)
                dist_creator.save(
                    os.path.join(
                        plot_path,
                        savefile_dist.format("real", source, detector, seed,
                                             ptype)))

    # plot assossciation skymap + corner plots, all sources, all particle types, only 2 seeds, TA only

    result_args_list = []
    for ptype in ptypes:
        for source in sources:
            for seed in seeds:

                fig_args = {
                    "model": "joint_gmf",
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

                output_figs = OutputFigures(fig_args)

                savefile_assos_skymap = os.path.join(
                    plot_path,
                    "tmp_GMFplots_{0}_{1}_{2}_{3}_{4}_{5}_assos_skymap.png".
                    format(*fig_args.values()))

                savefile_corner = os.path.join(
                    plot_path,
                    "tmp_GMFplots_{0}_{1}_{2}_{3}_{4}_{5}_corner.png".format(
                        *fig_args.values()))

                dominant_sources_dict = output_figs.association_skymap(
                    savefile_assos_skymap)

                output_figs.corner_data(savefile_corner)

                result_args = {
                    "model": "joint_gmf",
                    "dtype": "real",
                    "source": source,
                    "detector": detector,
                    "seed": seed,
                    "ptype": ptype,
                    "tight_B": False,
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