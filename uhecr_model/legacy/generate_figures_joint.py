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

    # for joint model only
    sources = ["SBG_23", "2FHL_250Mpc", "swift_BAT_213"]
    detectors = ["TA2015", "auger2014"]
    models = ["arrival_direction", "joint"]
    dtypes = ["real", "sim"]
    ptype = ["p", "He", "N", "O", "Si", "Fe"]
    seeds = [19990308, 16852056]

    output_file = "{0}_fit_{1}_{2}_{3}_{4}_{5}.h5"

    result_args_list = []

    for source in sources:
        for detector in detectors:
            for seed in seeds:
                # no need sim_model since for auger / TA results we didnt include them

                fig_args = {
                    "model": "joint",
                    "dtype": "real",
                    "source": source,
                    "detector": detector,
                    "seed": seed,
                    "ptype": "p",
                }
                print("config: " + "    ".join([
                    "{0} : {1}".format(key, value)
                    for key, value in fig_args.items()
                ]))

                output_figs = OutputFigures(fig_args)

                savefile_assos_skymap = os.path.join(
                    plot_path,
                    "TAplots_{0}_{1}_{2}_{3}_{4}_{5}_assos_skymap.png".format(
                        *fig_args.values()))

                savefile_corner = os.path.join(
                    plot_path,
                    "TAplots_{0}_{1}_{2}_{3}_{4}_{5}_corner.png".format(
                        *fig_args.values()))

                dominant_sources_dict = output_figs.association_skymap(
                    savefile_assos_skymap)

                output_figs.corner_data(savefile_corner)

                result_args = {
                    "model": "joint",
                    "dtype": "real",
                    "source": source,
                    "detector": detector,
                    "seed": seed,
                    "ptype": "p",
                    "tight_B": False,
                    "dominant_sources": list(dominant_sources_dict.keys()),
                    "counts": list(dominant_sources_dict.values())
                }

                result_args_list.append(result_args)

    # write dominant sources to some csv file
    field_names = [
        "model", "dtype", "source", "detector", "seed", "ptype", "tight_B",
        "dominant_sources", "counts"
    ]

    with open(os.path.join(result_path, "dominant_source_list.csv"), "w") as f:
        writer = csv.DictWriter(f, fieldnames=field_names)

        writer.writeheader()
        for result_args in result_args_list:
            writer.writerow(result_args)

    # TA, comparing models for each source and each seed (choose best seed)
    savefile_dist = "TAplots_{0}_{1}_{2}_{3}_{4}_dist.png"
    for source in sources:
        for seed in seeds:
            print("config: TA, all models, source: {0}, seed: {1}".format(
                source, seed))
            output_files = [
                os.path.join(
                    output_path,
                    output_file.format(model, "real", source, "TA2015", seed,
                                       "p")) for model in models
            ]
            source_label = "Swift-BAT" if source == "swift_BAT_213" else source.split(
                "_")[0]

            dist_creator = SourceUHECRDist()

            dist_creator.get_fs(output_files)

            dist_creator.plotdist(models, title=source_label)
            dist_creator.save(
                os.path.join(
                    plot_path,
                    savefile_dist.format("real", source, "TA2015", seed, "p")))

    # TA: TA vs Auger, all sources, all seeds
    savefile_dist = "TAplots_{0}_{1}_{2}_{3}_TAvsAuger_dist.png"
    detector_labels = ["TA", "Auger"]
    for source in sources:
        for seed in seeds:
            print("config: TA vs auger, source: {0}, seed: {1}".format(
                source, seed))
            output_files = [
                os.path.join(
                    output_path,
                    output_file.format("joint", "real", source, detector, seed,
                                       "p")) for detector in detectors
            ]
            source_label = "Swift-BAT" if source == "swift_BAT_213" else source.split(
                "_")[0]

            dist_creator = SourceUHECRDist()

            dist_creator.get_fs(output_files)

            dist_creator.plotdist(detector_labels, title=source_label)
            dist_creator.save(
                os.path.join(plot_path,
                             savefile_dist.format("real", source, seed, "p")))

    # TA: source comparison, all detectors, all seeds
    savefile_dist = "TAplots_{0}_{1}_{2}_{3}_allsrc_dist.png"
    src_labels = ["SBG", "2FHL", "Swift-BAT"]
    detector_labels = ["TA", "Auger"]
    for i, detector in enumerate(detectors):
        for seed in seeds:
            print("config: all sources, detector : {0}, seed: {1}".format(
                detector, seed))
            output_files = [
                os.path.join(
                    output_path,
                    output_file.format("joint", "real", source, detector, seed,
                                       "p")) for source in sources
            ]

            dist_creator = SourceUHECRDist()

            dist_creator.get_fs(output_files)

            dist_creator.plotdist(src_labels, title=detector_labels[i])
            dist_creator.save(
                os.path.join(plot_path,
                             savefile_dist.format("real", detector, seed,
                                                  "p")))
