'''
Python script to generate figures from output files.
Converted from figures.ipynb.
'''
import os
import argparse
import csv

from matplotlib.pyplot import plot

# from figures.output_figures import OutputFigures, SourceUHECRDist
from figures.deflection_figures import DeflectionFigures
from figures.gmfsim_figures import GMFSimFigures

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
    sources = ["SBG_23"]
    # detectors = ["TA2015", "auger2014"]
    models = ["arrival_direction", "joint", "joint_gmf"]
    # dtypes = ["real", "sim"]
    ptypes = ["p", "He", "N", "O", "Si", "Fe"]
    # seeds = [19990308, 16852056, 65492186, 9999999, 9953497]
    seeds = [19990308]

    output_file = "{0}_fit_{1}_{2}_{3}_{4}_{5}_tightB.h5"
    model_labels = ["Joint", "Joint + GMF"]

    # defl_fig_generator = DeflectionFigures()

    # nodefl_savefile = os.path.join(plot_path,
    #                                "GMFplots_nodefl_smear_skymap.png")

    # defl_fig_generator.plot_smearing(nodefl_savefile, ptype="p", nodefl=True)

    # for ptype in ptypes:
    #     print("Current ptype: ", ptype)
    #     defl_savefile = os.path.join(
    #         plot_path, "GMFplots_{0}_defl_smear_skymap.png".format(ptype))

    #     defl_fig_generator.plot_smearing(defl_savefile, ptype=ptype)

    #     # select index = 10, which is ~ near hotspot
    #     defl_selUHECR_savefile = os.path.join(
    #         plot_path,
    #         "GMFplots_{0}_selUHECR_defl_smear_skymap.png".format(ptype))
    #     defl_fig_generator.plot_smearing(defl_selUHECR_savefile,
    #                                      ptype=ptype,
    #                                      sel_uhecr_idx=10)

    #     # plot deflection circles
    #     circle_savefile = os.path.join(
    #         plot_path, "GMFplots_{0}_defl_circles_skymap.png".format(ptype))

    #     defl_fig_generator.plot_circles(circle_savefile, ptype=ptype)

    # plot GMF sim plots
    sim_output_file = "{0}_sim_{1}_{2}_{3}_{4}_{5}.h5"
    savefile = "GMFplots_{0}_{1}_{2}_{3}_{4}_{6}_{5}sim.png"

    for source in sources:
        for seed in seeds:
            for ptype in ptypes:
                print("GMF sims: {0}, {1}, {2}".format(source, seed, ptype))

                sim_output_file = os.path.join(
                    output_path,
                    sim_output_file.format("joint_gmf", source, "TA2015", seed,
                                           ptype, "notightB"))

                gmfplots_creator = GMFSimFigures(sim_output_file)

                # create uhecr skymap of randomized particles at galactic boundary
                gmfplots_creator.uhecr_gb_skymap(
                    os.path.join(
                        plot_path,
                        savefile.format("joint_gmf", source, "TA2015", seed,
                                        ptype, "notightB", "uhecr_gb_skymap")))

                # unlensed and lensed maps
                gmfplots_creator.unlensed_skymap(
                    os.path.join(
                        plot_path,
                        savefile.format("joint_gmf", source, "TA2015", seed,
                                        ptype, "notightB", "unlensed_skymap")))

                gmfplots_creator.lensed_skymap(
                    os.path.join(
                        plot_path,
                        savefile.format("joint_gmf", source, "TA2015", seed,
                                        ptype, "notightB", "lensed_skymap")))

                # simulated UHECR
                gmfplots_creator.sim_uhecr_skymap(
                    os.path.join(
                        plot_path,
                        savefile.format("joint_gmf", source, "TA2015", seed,
                                        ptype, "notightB",
                                        "sim_uhecr_skymap")))
