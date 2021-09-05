'''Collection of functions that generate figures from the fit outputs. Creates figures from fit results using both simulated and real datasets.'''

from fancy.analysis import results
import numpy as np
import os
import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns
from pandas import DataFrame
from astropy.coordinates import SkyCoord
from astropy import units as u
import csv
from copy import deepcopy

from fancy import Data, Results
from fancy.plotting import Corner
from fancy.plotting.allskymap_cartopy import AllSkyMapCartopy as AllSkyMap
from fancy.plotting.colours import *
from fancy.interfaces.stan import Direction

# use minimalist style
plt.style.use("minimalist")

# paths to important files
path_to_this_file = os.path.abspath(os.path.dirname(__file__))
stan_path = os.path.join(path_to_this_file, "..", "stan")
source_file = os.path.join(path_to_this_file, "..", "data", "sourcedata.h5")
uhecr_file = os.path.join(path_to_this_file, "..", "data", "UHECRdata.h5")
table_path = os.path.join(path_to_this_file, "..", "tables")
output_path = os.path.join(path_to_this_file, "..", "outputs")
output_sim_path = os.path.join(output_path, "sim")
output_data_path = os.path.join(output_path, "data")

plot_path = os.path.join(path_to_this_file, "..", "..", "..", "plots")

# make corresponding directories if they dont exist
if not os.path.exists(plot_path):
    os.mkdir(plot_path)


class OutputFigures():
    def __init__(self, fig_args, sim_model=None, header="tmp", verbose=False):
        '''Collection of figures obtained from output of fitting.'''
        # parameters set from dictionary
        self.source = fig_args["source"]
        self.detector = fig_args["detector"]
        self.ptype = fig_args["ptype"]
        self.model = fig_args["model"]

        self.seed = fig_args["seed"]
        self.end_label = fig_args["end_label"]

        # store the dictionary too (for dominant sources)
        self.fig_args = fig_args

        # simulation model and choice if we want to plot simulation or not
        self.sim = False if sim_model is None else True  # enable simulation or not
        self.sim_model = sim_model

        self.verbose = verbose

        self._initialize_files(header)

        # obtain detector properties / params from imports
        self.detector_properties, self.detector_params = self._get_detectorimports(
        )

        # initialize the data object and add relevant data into it
        self._initialize_data()

        # format style for skymap labels (from "mpl" or "TA")
        self.skymap_label_fmt = "default"

    def _initialize_files(self, header):
        '''Initialize filenames'''
        # # output file from simulation / fitting

        if self.end_label is not None:
            output_file = "{0}_{1}_{2}_{3}_{4}.h5".format(
                self.source, self.detector, self.seed, self.ptype,
                self.end_label)
            # create savefile template
            savefile_templ = "{0}_{1}_{2}_{3}_{4}_{{0}}.png".format(
                self.source, self.detector, self.seed, self.ptype,
                self.end_label)
        else:
            output_file = "{0}_{1}_{2}_{3}.h5".format(self.source,
                                                      self.detector, self.seed,
                                                      self.ptype)
            savefile_templ = "{0}_{1}_{2}_{3}_{{0}}.png".format(
                self.source, self.detector, self.seed, self.ptype)

        # create header directories
        header_path = os.path.join(plot_path, header)
        if not os.path.exists(header_path):
            os.mkdir(header_path)

        # create subdirectories that contain simulation / data results
        for ftype in ["simulation", "data"]:
            if not os.path.exists(os.path.join(header_path, ftype)):
                os.mkdir(os.path.join(header_path, ftype))

            if self.sim:  # then create subdirectories of simulation models
                if not os.path.exists(
                        os.path.join(header_path, ftype, self.sim_model)):
                    os.mkdir(os.path.join(header_path, ftype, self.sim_model))

                if not os.path.exists(
                        os.path.join(header_path, ftype, self.sim_model,
                                     self.model)):
                    os.mkdir(
                        os.path.join(header_path, ftype, self.sim_model,
                                     self.model))

            else:  # for data, no need this extra layer
                if not os.path.exists(
                        os.path.join(header_path, ftype, self.model)):
                    os.mkdir(os.path.join(header_path, ftype, self.model))

        if self.sim:

            self.sim_output_file = os.path.join(output_sim_path,
                                                self.sim_model, "simulation",
                                                output_file)

            self.analysis_output_file = os.path.join(output_sim_path,
                                                     self.sim_model, "fit",
                                                     self.model, output_file)

            self.savefile_templ = os.path.join(header_path, "simulation",
                                               self.sim_model, self.model,
                                               savefile_templ)

        else:
            self.analysis_output_file = os.path.join(output_data_path,
                                                     self.model, output_file)

            self.savefile_templ = os.path.join(header_path, "data", self.model,
                                               savefile_templ)

        # where to write dominant sources to, some csv file
        self.dominant_sources_file = os.path.join(header_path,
                                                  "dominant_sources.csv")

        if self.verbose:
            print("Analysis Output File: ", self.analysis_output_file)
            print("Saving to: ", self.savefile_templ)

            if self.sim:
                print("Reading from Simulation Output File: ",
                      self.sim_output_file)

    def _get_detectorimports(self):
        '''Get variables imported by (detector_name).py'''
        if self.detector == "TA2015":
            from fancy.detector.TA2015 import detector_properties, detector_params
        elif self.detector == "auger2014":
            from fancy.detector.auger2014 import detector_properties, detector_params
        elif self.detector == "auger2010":
            from fancy.detector.auger2010 import detector_properties, detector_params
        else:
            raise Exception("Undefined detector type!")

        return detector_properties, detector_params

    def _initialize_data(self):
        '''Initialize data object based on sim / real data'''
        self.data = Data()

        if self.sim:
            self.data.from_file(self.sim_output_file)

        else:
            self.data.add_source(source_file, self.source)
            self.data.add_uhecr(uhecr_file, self.detector)
            self.data.add_detector(self.detector_properties)

    def src_uhecr_skymap(self, coord="G", exposure="map"):
        '''
        Plot skymap with sources + UHECR, obtained from Data object.
        Optionally add either exposure map or exposure limit.
        '''
        # read in initializations
        omega_src = Direction(self.data.source.unit_vector)
        omega_arr = Direction(self.data.uhecr.unit_vector)
        energy_arr = self.data.uhecr.energy

        # convert source / arrival directions to equatorial / galactic
        if coord == "G":
            x_src, y_src = 180 - omega_src.glons, omega_src.glats
            x_arr, y_arr = 180 - omega_arr.glons, omega_arr.glats
        elif coord == "E":
            x_src, y_src = 180 - omega_src.ras, omega_src.decs
            x_arr, y_arr = 180 - omega_arr.ras, omega_arr.decs

        Eth = self.data.detector.Eth
        Emax = np.ceil(np.max(energy_arr) / 10.) * 10.

        uhecr_color = [lightblue, midblue, darkblue]
        uhecr_cmap = mpl.colors.ListedColormap(uhecr_color)
        energy_bins = np.logspace(np.log(Eth), np.log(Emax), 4, base=np.e)
        uhecr_norm = mpl.colors.BoundaryNorm(energy_bins, uhecr_cmap.N)

        # Legend
        legend_elements = [
            mpl.lines.Line2D([0], [0],
                             marker='o',
                             color='w',
                             label='sources',
                             markersize=10,
                             markerfacecolor='k'),
            mpl.lines.Line2D([0], [0],
                             marker='o',
                             color='w',
                             label='UHECRs',
                             markersize=15,
                             markerfacecolor=midblue,
                             alpha=0.8)
        ]

        # create skymap
        skymap = AllSkyMap(projection='moll', lon_0=180)
        skymap.set_gridlines(label_fmt=self.skymap_label_fmt)

        # sources
        skymap.scatter(x_src, y_src, s=10.0, color='k', alpha=1.0, zorder=5)

        # UHECRs
        for lon, lat, E in np.nditer([x_arr, y_arr, energy_arr]):
            i = np.digitize(E, energy_bins) - 1
            skymap.tissot(lon,
                          lat,
                          3.0 + (i * 2),
                          30,
                          facecolor=uhecr_cmap.colors[i],
                          alpha=0.8,
                          zorder=i + 2)

        # exposure
        if exposure == "map":
            skymap.exposure_map(self.data.detector.params, coord=coord)
        elif exposure == "limit":
            skymap.exposure_limit(self.data.detector.limiting_dec.deg,
                                  coord=coord,
                                  s=2,
                                  color=grey,
                                  alpha=1)
        else:
            raise ValueError(
                "Exposure plot type {0} not defined.".format(exposure))

        # Annotations
        skymap.legend(handles=legend_elements,
                      loc='upper right',
                      bbox_to_anchor=(1., 1.),
                      fontsize=16,
                      fancybox=True)
        # skymap.title("{0} + {1}".format(src_label, uhecr_label))

        # Colorbar
        cb_ax = plt.axes([0.25, 0.07, .5, .05], frameon=False)
        bar = mpl.colorbar.ColorbarBase(cb_ax,
                                        norm=uhecr_norm,
                                        cmap=uhecr_cmap,
                                        orientation='horizontal',
                                        drawedges=True,
                                        alpha=1)
        bar.set_label('$\hat{E}$ / EeV', color='k', fontsize=16)
        bar.ax.tick_params(labelsize=16)

        skymap.save(self.savefile_templ.format("skymap"))

    def eval_association_probs(self):
        '''
        Evaluate the assocation probabilities between sources and UHECR
        from some output file
        '''
        # TODO: check how lp is derived with arrival direction model in analysis / stan
        # TODO: add some way to write to dominant sources list
        #       - simply append to a given .csv file in plot_path/header, written in
        # Log probability
        results = Results(self.analysis_output_file)
        keys = ['lp']
        chain = results.get_chain(keys)
        logprob = chain['lp'].transpose(1, 2, 0)
        N = np.shape(logprob)[0]

        # Account for background component
        Ns = np.shape(logprob)[1] - 1

        # Calculate association probabilities for each source-UHECR combo
        uhecr_p = []
        for lp in logprob:
            lps = []
            for src in range(Ns + 1):
                lps.append(np.mean(np.exp(lp[src])))

            norm = sum(lps)
            ps = []
            for src in range(Ns + 1):
                ps.append(lps[src] / norm)
            uhecr_p.append(ps)

        # Normalise line weights
        pmax = max(max(uhecr_p))

        # Find names of dominant sources
        threshold_probability = 0.1

        dominant = []
        for p in uhecr_p:
            # for i in range(data.source.N):
            for i in range(Ns):
                if p[i] > threshold_probability:
                    dominant.append(i)

        seen = set()
        inds = []
        for d in dominant:
            if d not in seen:
                inds.append(d)
                seen.add(d)

        # dominant_sources = [self.data.source.name[i] for i in inds]

        if not self.sim:
            # sort so that those with largest assossations appear first
            N_assos = {}
            for i in inds:
                N_assos[self.data.source.name[i].decode("UTF-8")] = len(
                    np.argwhere([d == i for d in dominant]))

            N_assos_sorted_dict = {
                k: v
                for k, v in sorted(N_assos.items(), key=lambda item: item[1])
                [::-1]
            }

            # write to csv file
            self._write_dom_srcs(N_assos_sorted_dict)

        # print("Dominant sources: ", [self.data.source.name[i] for i in inds])

        return uhecr_p

    def _write_dom_srcs(self, dom_srcs_dict):
        '''Write dominant sources to the csv files'''

        fig_args_dom_srcs = deepcopy(self.fig_args)

        fig_args_dom_srcs.update({
            "dominant_sources": list(dom_srcs_dict.keys()),
            "counts": list(dom_srcs_dict.values()),
            "simulation": self.sim
        })

        field_names = [
            "model", "simulation", "source", "detector", "seed", "ptype",
            "end_label", "dominant_sources", "counts"
        ]

        # flag so that we dont write the header files all the time.
        file_exists = os.path.exists(self.dominant_sources_file)

        with open(self.dominant_sources_file, "a") as f:
            writer = csv.DictWriter(f, fieldnames=field_names)

            # only write header if the file has just been created
            if not file_exists:
                writer.writeheader()
            writer.writerow(fig_args_dom_srcs)

    def association_skymap(self, coord="G"):
        '''Plot association skymap between sources and UHECRs'''

        uhecr_p = self.eval_association_probs()

        # labels
        src_label = self.data.source.label.decode("UTF-8") if isinstance(
            self.data.source.label, bytes) else self.data.source.label

        # print(type(self.data.source.label))

        src_label = "Swift-BAT" if src_label == "swift_BAT_213" else src_label.split(
            "_")[0]

        # uhecr_label = self.data.uhecr.label
        # read in initializations
        omega_src = Direction(self.data.source.unit_vector)
        omega_arr = Direction(self.data.uhecr.unit_vector)
        energy_arr = self.data.uhecr.energy

        # convert source / arrival directions to equatorial / galactic
        if coord == "G":
            x_src, y_src = 180. - omega_src.glons, omega_src.glats
            x_arr, y_arr = 180. - omega_arr.glons, omega_arr.glats
        elif coord == "E":
            x_src, y_src = omega_src.ras, omega_src.decs
            x_arr, y_arr = omega_arr.ras, omega_arr.decs

        Eth = self.data.detector.Eth
        Emax = np.ceil(np.max(energy_arr) / 10.) * 10.

        uhecr_color = [lightblue, midblue, darkblue]
        uhecr_cmap = mpl.colors.ListedColormap(uhecr_color)
        energy_bins = np.logspace(np.log(Eth), np.log(Emax), 4, base=np.e)
        uhecr_norm = mpl.colors.BoundaryNorm(energy_bins, uhecr_cmap.N)

        # Legend
        legend_elements = [
            mpl.lines.Line2D([0], [0],
                             marker='o',
                             color='w',
                             label=src_label,
                             markersize=10,
                             markerfacecolor='k'),
            mpl.lines.Line2D([0], [0],
                             marker='o',
                             color='w',
                             label='UHECRs',
                             markersize=15,
                             markerfacecolor=midblue,
                             alpha=0.8)
        ]

        # plot
        skymap = AllSkyMap(projection='moll', lon_0=180)
        skymap.set_gridlines(label_fmt=self.skymap_label_fmt)

        # sources
        skymap.scatter(x_src, y_src, s=10.0, color='k', alpha=1.0, zorder=5)

        # UHECRs
        for lon, lat, E in np.nditer([x_arr, y_arr, energy_arr]):
            i = np.digitize(E, energy_bins) - 1
            skymap.tissot(lon,
                          lat,
                          3.0 + (i * 2),
                          30,
                          facecolor=uhecr_cmap.colors[i],
                          alpha=0.8,
                          zorder=i + 2)

        # Association
        # add some way to include dominant sources with this association plot later
        pmax = np.max(np.max(uhecr_p))

        for i, p in enumerate(uhecr_p):
            for j, psrc in enumerate(p[0:self.data.source.N]):
                if psrc > 0.001:
                    skymap.geodesic(x_arr[i],
                                    y_arr[i],
                                    x_src[j],
                                    y_src[j],
                                    color='k',
                                    lw=3,
                                    alpha=psrc / pmax,
                                    zorder=10)

        # Annotations
        skymap.legend(handles=legend_elements,
                      loc='upper right',
                      bbox_to_anchor=(1.1, 1.),
                      fontsize=16,
                      fancybox=True)

        # exposure limit
        skymap.exposure_limit(
            self.data.detector.limiting_dec.deg,
            coord=coord,
            s=0.5,
            marker="o",
            color=lightpurple,
            alpha=0.01,
        )

        # Colorbar
        cb_ax = plt.axes([0.25, 0.07, .5, .05], frameon=True)
        bar = mpl.colorbar.ColorbarBase(cb_ax,
                                        norm=uhecr_norm,
                                        cmap=uhecr_cmap,
                                        orientation='horizontal',
                                        drawedges=True,
                                        alpha=1)
        bar.set_label('$\hat{E}$ / EeV', color='k', fontsize=16)
        bar.ax.tick_params(labelsize=16)

        skymap.save(self.savefile_templ.format("assos_skymap"))

    def corner(self):
        '''Plot corner plot'''
        if self.sim:
            self.corner_sim()
        else:
            self.corner_data()

    def corner_data(self):
        '''Plot corner plot of data'''
        # Get chains from joint fit
        results_fit = Results(self.analysis_output_file)

        # get keys and corresponding labels
        labels = {}
        if self.model.find("join") != -1:
            keys = ['alpha', 'B', 'f', 'L']

            labels['B'] = r'$B$ / $\mathrm{nG}$'
            labels['alpha'] = r'$\alpha$'
            labels['f'] = r'$f$'
            labels['L'] = r'$L$ / $10^{38}$ $\mathrm{yr}^{-1}$'

            chain = results_fit.get_chain(keys)
            chain["L"] = chain["L"] * 10

        elif self.model == "arrival_direction":
            keys = ['kappa', 'L', 'f']

            labels['L'] = r'$L$ / $10^{38}$ $\mathrm{yr}^{-1}$'
            labels['kappa'] = r'$\kappa$'
            labels['f'] = r'$f$'

            chain = results_fit.get_chain(keys)
            # for better scaling
            chain["L"] = chain["L"] * 10  # 10^-38 yr^-1
            chain["kappa"] = np.log10(chain["kappa"])

        # Make nicely labelled dict
        chain_for_df = {}
        for key in keys:
            chain_for_df[labels[key]] = chain[key]

        # Make ordered dataframe
        df = DataFrame(data=chain_for_df)
        df = df[[labels[key] for key in keys]]

        corner = Corner(df,
                        color=midblue,
                        contour_color=midblue_contour,
                        end_label=self.end_label)
        corner.save(self.savefile_templ.format("corner"))

    def corner_sim(self, cumul=False, seeds=None):

        if cumul:
            self.corner_sim_cumul(seeds)
        else:
            results_fit = Results(self.analysis_output_file)

            results_sim = Results(self.sim_output_file)
            truth_keys = ['F0', 'L', 'alpha', 'B', 'f']
            truth = results_sim.get_truths(truth_keys)
            info_keys = ['Eth', 'Eth_sim']
            info = results_sim.get_truths(info_keys)

            # Correct for different Eth in sim and fit
            # Also scale to plot units
            flux_scale = (info['Eth'] / info['Eth_sim'])**(1 - truth['alpha'])
            truth['F0'] = truth['F0'] * flux_scale  # km^-2 yr^-1
            truth[
                'L'] = truth['L'][0] * flux_scale / 1.0e39 * 10  # 10^-38 yr^-1

            if self.model == "arrival_direction":
                self.corner_sim_arrival(results_fit, truth)
            else:
                self.corner_sim_joint(results_fit, truth)

    def corner_sim_arrival(self, results_fit, truth):
        '''Corner plot for Arrival Direction model'''

        keys = ['kappa', 'L', 'f']
        chain = results_fit.get_chain(keys)

        chain['L'] = chain['L'] * 10  # 10^-38 yr^-1

        labels = {}
        labels['L'] = r'$L$ / $10^{38}$ $\mathrm{yr}^{-1}$'
        labels['kappa'] = r'$\kappa$'
        labels['f'] = r'$f$'

        truths = [truth["L"], truth["f"]]

        # Make nicely labelled dict
        chain_for_df = {}
        for key in keys:
            chain_for_df[labels[key]] = chain[key]

        # Make ordered dataframe
        df = DataFrame(data=chain_for_df)
        df = df[[labels['L'], labels['f'], labels['kappa']]]

        corner = Corner(df,
                        truths,
                        color=midblue,
                        contour_color=midblue_contour,
                        end_label=self.end_label)

        corner.save(self.savefile_templ.format("corner"))

    def corner_sim_joint(self, results_fit, truth):
        '''Plot corner plot resulting from simulations'''

        # keys = ['F0', 'L', 'alpha', 'B', 'f']
        keys = ['L', 'alpha', 'B', 'f']
        chain = results_fit.get_chain(keys)

        # Convert form Stan units to plot units
        # chain['F0'] = chain['F0'] / 1.0e3  # km^-2 yr^-1
        chain['L'] = chain['L'] * 10  # 10^-38 yr^-1

        labels = {}
        labels['L'] = r'$L$ / $10^{38}$ $\mathrm{yr}^{-1}$'
        # labels['F0'] = r'$F_0$ / $\mathrm{km}^{-2} \ \mathrm{yr}^{-1}$'
        labels['B'] = r'$B$ / $\mathrm{nG}$'
        labels['alpha'] = r'$\alpha$'
        labels['f'] = r'$f$'

        params = np.column_stack([chain[key] for key in keys])
        truths = [truth[key] for key in keys]

        # Make nicely labelled dict
        chain_for_df = {}
        for key in keys:
            chain_for_df[labels[key]] = chain[key]

        # Make ordered dataframe
        df = DataFrame(data=chain_for_df)
        # df = df[[
        #     labels['F0'], labels['L'], labels['alpha'], labels['B'],
        #     labels['f']
        # ]]
        df = df[[labels['L'], labels['alpha'], labels['B'], labels['f']]]

        corner = Corner(df,
                        truths,
                        color=purple,
                        contour_color=purple_contour,
                        end_label=self.end_label)
        corner.save(self.savefile_templ.format("corner"))

    def _get_files_for_corner(self, seed):
        '''Get output files used for averaged out corner plot'''
        if self.end_label is not None:
            analysis_output_file = "{0}_fit_{1}_{2}_{3}_{4}_{5}.h5".format(
                self.model, self.source, self.detector, seed, self.ptype,
                self.end_label)
        else:
            analysis_output_file = "{0}_fit_{1}_{2}_{3}_{4}.h5".format(
                self.model, self.source, self.detector, seed, self.ptype)

        return analysis_output_file

    def corner_sim_cumul(self, seeds):
        '''Cumulative distribution of corner plots from simulation'''

        # keys = ['F0', 'L', 'alpha', 'B', 'f']
        keys = ['L', 'alpha', 'B', 'f']
        chain_avgs = {key: 0 for key in keys}
        chain_list = []
        Nseeds = len(seeds)

        # get chains for each output file
        # ignores self.outpuf_file
        for seed in seeds:
            output_file = self._get_files_for_corner(seed)
            chain = Results(output_file).get_chain(keys)

            chain_list.append(chain)

        # evaluate averaged value for each seed
        for key in keys:
            chain_sum = 0
            for i in range(Nseeds):
                chain_sum += chain_list[i][key]

            chain_sum /= Nseeds

            chain_avgs[key] = chain_sum

        # Convert form Stan units to plot units
        # chain_avgs['F0'] = chain_avgs['F0'] / 1.0e3  # km^-2 yr^-1
        chain_avgs['L'] = chain_avgs['L'] * 10  # 10^-38 yr^-1

        # Get truths from simulation
        results_sim = Results(self.sim_output_file)

        # truth_keys = ['F0', 'L', 'alpha', 'B', 'f']
        truth_keys = keys
        truth = results_sim.get_truths(truth_keys)
        info_keys = ['Eth', 'Eth_sim']
        info = results_sim.get_truths(info_keys)

        # Correct for different Eth in sim and fit
        # Also scale to plot units
        flux_scale = (info['Eth'] / info['Eth_sim'])**(1 - truth['alpha'])
        # truth['F0'] = truth['F0'] * flux_scale  # km^-2 yr^-1
        truth['L'] = truth['L'][0] * flux_scale / 1.0e39 * 10  # 10^-38 yr^-1

        labels = {}
        labels['L'] = r'$L$ / $10^{38}$ $\mathrm{yr}^{-1}$'
        # labels['F0'] = r'$F_0$ / $\mathrm{km}^{-2} \ \mathrm{yr}^{-1}$'
        labels['B'] = r'$B$ / $\mathrm{nG}$'
        labels['alpha'] = r'$\alpha$'
        labels['f'] = r'$f$'

        params = np.column_stack([chain_avgs[key] for key in keys])
        truths = [truth[key] for key in keys]

        # Make nicely labelled dict
        chain_for_df = {}
        for key in keys:
            chain_for_df[labels[key]] = chain_avgs[key]

        # Make ordered dataframe
        df = DataFrame(data=chain_for_df)
        # df = df[[
        #     labels['F0'], labels['L'], labels['alpha'], labels['B'],
        #     labels['f']
        # ]]
        df = df[[labels['L'], labels['alpha'], labels['B'], labels['f']]]

        corner = Corner(df, truths, color=purple, contour_color=purple_contour)
        corner.save(self.savefile_templ.format("corner_cumul"))

    def dist(self, param="f"):
        '''1-D Distribution of param similar to f-distribution.'''

        fig, ax = plt.subplots(figsize=(6, 4))

        fig_params = {
            "f": ["$f$", [0.9, 1], [0.9, 0.95, 1], "f"],
            "B": ["$B$ / nG", [0, 50], [0, 25, 50], "B", [0, 0.05, 0.1]],
            "alpha": ["$\alpha$", [3, 6], [3, 4.5, 6], "\alpha"],
            "L": [r'$L$ / $10^{38}$ $\mathrm{yr}^{-1}$', 0, 0, "L"]
        }

        if self.end_label == "tight_B" or self.end_label == "limitL":
            fig_params["B"] = [
                "$B$ / nG", [0, 10], [0, 5.0, 10], "B",
                [0, 0.25, 0.5, 0.75, 1.]
            ]

        result_fit = Results(self.analysis_output_file).get_chain([param
                                                                   ])[param]

        if param == "L":
            result_fit = result_fit * 10

        sns.distplot(result_fit,
                     hist=False,
                     kde_kws={
                         'shade': True,
                         'lw': 2
                     },
                     color=purple)

        if self.sim:
            f_true = Results(self.sim_output_file).get_truths(['f'])['f']
            ax.axvline(f_true, 0, 10, color='k', zorder=3, lw=2., alpha=0.7)

        # ax.set_xlim(fig_params[param][1])
        ax.set_xlabel(fig_params[param][0], fontsize=14)
        # ax.set_xticks(fig_params[param][2])
        # ax.set_yticks(fig_params[param][4])
        ax.set_ylabel('$P({0} | \hat{{E}}, \hat{{\omega}})$'.format(
            fig_params[param][3]),
                      fontsize=14)
        ax.tick_params(axis="x", labelsize=14)
        ax.tick_params(axis="y", labelsize=14)

        fig.savefig(self.savefile_templ.format("{0}_dist".format(param)),
                    bbox_inches='tight',
                    dpi=300)


class SourceUHECRDist():
    def __init__(self,
                 fig_args,
                 config_dict,
                 sim_model=None,
                 header="tmp",
                 verbose=False):
        '''Class that organizes plotting of source-UHECR association fraction distribution'''

        self.source = fig_args["source"]
        self.detector = fig_args["detector"]
        self.ptype = fig_args["ptype"]
        self.model = fig_args["model"]

        self.seed = fig_args["seed"]
        self.end_label = fig_args["end_label"]

        # container that contains the configuration used to plot against
        self.config_dict = config_dict

        self.sim = False if sim_model is None else True  # enable simulation or not
        self.sim_model = sim_model

        self.verbose = verbose
        self.header = header

        self.color_list = [grey, lightpurple, lightblue]

    def _get_output_fmt_list(self, compare_label, value):
        '''Get list of values used for string format in output file'''
        if compare_label == "model":
            return [self.source, self.detector, self.seed, self.ptype]
        elif compare_label == "source":
            return [value, self.detector, self.seed, self.ptype]
        elif compare_label == "detector":
            return [self.source, value, self.seed, self.ptype]
        elif compare_label == "seed":
            return [self.source, self.detector, value, self.ptype]
        elif compare_label == "ptype":
            return [self.source, self.detector, self.seed, value]

    def _get_savefile_fmt_list(self, compare_label):
        '''Get list of values used for string format in output file'''
        if compare_label == "model":
            return [
                self.source, self.detector, self.seed, self.ptype,
                compare_label
            ]
        elif compare_label == "source":
            return [
                self.model, self.detector, self.seed, self.ptype, compare_label
            ]
        elif compare_label == "detector":
            return [
                self.model, self.source, self.seed, self.ptype, compare_label
            ]
        elif compare_label == "seed":
            return [
                self.model, self.source, self.detector, self.ptype,
                compare_label
            ]
        elif compare_label == "ptype":
            return [
                self.model, self.source, self.detector, self.seed,
                compare_label
            ]

    def _create_savefile_templ(self, compare_label, cumul):
        '''Get save file template'''
        # str list for savefile formatting
        # contains everything without the ones we are comparing with
        savefile_fmt_list = self._get_savefile_fmt_list(compare_label)
        # savefile_fmt_list.append(compare_label)
        if cumul:
            if compare_label == "ptype":
                del savefile_fmt_list[3]
            elif compare_label == "seed":
                raise Exception("Cannot used seed with cumul distributioN!")
            else:
                del savefile_fmt_list[2]

            if self.end_label is None:
                savefile_templ = "{0}_{1}_{2}_all{3}dist_avg.png".format(
                    *savefile_fmt_list)
            else:
                savefile_fmt_list.append(self.end_label)
                savefile_templ = "{0}_{1}_{2}_{4}_all{3}dist_avg.png".format(
                    *savefile_fmt_list)
        else:
            if self.end_label is None:
                savefile_templ = "{0}_{1}_{2}_{3}_all{4}dist.png".format(
                    *savefile_fmt_list)
            else:
                savefile_fmt_list.append(self.end_label)
                savefile_templ = "{0}_{1}_{2}_{3}_{5}_all{4}dist.png".format(
                    *savefile_fmt_list)

        return savefile_templ

    def _initialize_savefiles(self, compare_label, cumul=False):
        '''Initialize savefiles'''

        # create savefile template
        savefile_templ = self._create_savefile_templ(compare_label, cumul)

        # below is new code.
        header_path = os.path.join(plot_path, self.header)
        if not os.path.exists(header_path):
            os.mkdir(header_path)

        # create subdirectories that contain simulation / data results
        for ftype in ["simulation", "data"]:
            if not os.path.exists(os.path.join(header_path, ftype)):
                os.mkdir(os.path.join(header_path, ftype))

            if self.sim:

                if not os.path.exists(
                        os.path.join(header_path, ftype, "fdists")):
                    os.mkdir(os.path.join(header_path, ftype, "fdists"))

                if not os.path.exists(
                        os.path.join(header_path, ftype, "fdists",
                                     self.sim_model)):
                    os.mkdir(
                        os.path.join(header_path, ftype, "fdists",
                                     self.sim_model))
            else:
                if not os.path.exists(
                        os.path.join(header_path, ftype, "fdists")):
                    os.mkdir(os.path.join(header_path, ftype, "fdists"))
        # create savefile
        if self.sim:
            self.savefile_templ = os.path.join(header_path, "simulation",
                                               "fdists", self.sim_model,
                                               savefile_templ)
        else:
            self.savefile_templ = os.path.join(header_path, "data", "fdists",
                                               savefile_templ)

    def _initialize_output_files(self, compare_label):
        '''Initialize output files '''

        # return output_files
        # get the batch of output file to read the source association fractions from
        output_files = []
        # if self.end_label is not None:
        #     output_file = "{0}_{1}_{2}_{3}_{4}.h5".format(
        #         self.source, self.detector, self.seed, self.ptype,
        #         self.end_label)
        # else:
        #     output_file = "{0}_{1}_{2}_{3}.h5".format(self.source,
        #                                               self.detector, self.seed,
        #                                               self.ptype)

        if self.end_label is not None:
            output_file = "{{0}}_{{1}}_{{2}}_{{3}}_{0}.h5".format(
                self.end_label)
        else:
            output_file = "{0}_{1}_{2}_{3}.h5"

        for value in list(self.config_dict[compare_label]):
            # str list for output file formatting
            # contains everything including the ones we are comparing with
            output_fmt_list = self._get_output_fmt_list(compare_label, value)

            # output file from simulation / fitting
            analysis_output_file = output_file.format(*output_fmt_list)

            model_name = self.model if compare_label != "model" else value

            if self.sim:
                # sim_fmt_list = output_fmt_list[
                #     1:]  # model information contained with dir structure
                sim_fmt_list = output_fmt_list
                sim_output_file = output_file.format(*sim_fmt_list)

                self.sim_output_file = os.path.join(output_sim_path,
                                                    self.sim_model,
                                                    "simulation",
                                                    sim_output_file)

                analysis_output_file = os.path.join(output_sim_path,
                                                    self.sim_model, "fit",
                                                    model_name,
                                                    analysis_output_file)

            else:
                analysis_output_file = os.path.join(output_data_path,
                                                    model_name,
                                                    analysis_output_file)

            output_files.append(analysis_output_file)
        return output_files

    def _initialize_files(self, compare_label, cumul=False):
        '''Initialize filenames'''

        self._initialize_savefiles(compare_label, cumul)

        if cumul:
            output_files = self._get_cumul_output_files(compare_label)
        else:
            output_files = self._initialize_output_files(compare_label)

        if self.verbose:
            print("Analysis Output Files: ", output_files)
            print("Saving to: ", self.savefile_templ)

            if self.sim:
                print("Reading from Simulation Output File: ",
                      self.sim_output_file)

        return output_files

    def _get_cumul_output_files(self, compare_label):
        '''Get output files for averaged distribution'''
        output_files = []
        if self.end_label is not None:
            output_file = "{{0}}_{{1}}_{{2}}_{{3}}_{0}.h5".format(
                self.end_label)
        else:
            output_file = "{0}_{1}_{2}_{3}.h5"

        for value in list(self.config_dict[compare_label]):
            output_files_per_value = []
            for seed in list(self.config_dict["seed"]):

                # str list for output file formatting
                # contains everything including the ones we are comparing with
                output_fmt_list = self._get_output_fmt_list(
                    compare_label, value)

                # modify seed part (index 3) to the seed in config dict
                output_fmt_list[2] = seed
                # print(output_fmt_list)
                # output_fmt_list_per_seed = deepcopy(output_fmt_list)
                # output_fmt_list_per_seed[3] = seed

                # # output file from simulation / fitting
                analysis_output_file = output_file.format(*output_fmt_list)

                # print(analysis_output_file)

                model_name = self.model if compare_label != "model" else value

                if self.sim:
                    sim_fmt_list = output_fmt_list  # model information contained with dir structure
                    sim_output_file = output_file.format(*sim_fmt_list)

                    self.sim_output_file = os.path.join(
                        output_sim_path, self.sim_model, "simulation",
                        sim_output_file)

                    analysis_output_file = os.path.join(
                        output_sim_path, self.sim_model, "fit", model_name,
                        analysis_output_file)

                output_files_per_value.append(analysis_output_file)

            # print(output_files_per_value)

            output_files.append(output_files_per_value)

            # print(output_files)

        return output_files

    def get_fs(self, compare_label, cumul=False):
        '''Get source-UHECR association fraction from output files'''

        output_files = self._initialize_files(compare_label, cumul=cumul)

        # print(output_files)

        self.f_list = []

        if cumul:  # evaluate cumulative distribution
            # then each element in fname_list must be a list of fnames with different models
            # each of this list will contain fnames with different seeds

            for output_files_per_model in output_files:
                f_list = []
                for output_file in output_files_per_model:
                    f_i = Results(output_file).get_chain(['f'])['f']
                    f_list.append(f_i)
                # print(np.array(f_list).shape)

                f_avg = np.mean(np.array(f_list), axis=0)
                print(f_avg.shape)
                self.f_list.append(f_avg)

        else:
            for output_file in output_files:
                # print(output_file)
                f_i = Results(output_file).get_chain(['f'])['f']
                self.f_list.append(f_i)

        return self.f_list

    def _get_labels(self, compare_label):
        '''Get labels used for plotting. For seed and ptype, use the list given'''
        if compare_label == "model":
            return ["Arrival Direction", "Joint", "Joint + GMF"]
        elif compare_label == "source":
            return ["SBG", "2FHL", "Swift-BAT"]
        elif compare_label == "detector":
            return ["TA", "Auger"]
        elif compare_label == "seed" or compare_label == "ptype":
            return None

    def plotdist(self,
                 compare_label="model",
                 cumul=False,
                 labels=None,
                 extend=False,
                 title=None):
        '''Plot the distribution'''

        fig, ax = plt.subplots(figsize=(6, 4))

        # get source fractions for the given output files
        self.get_fs(compare_label=compare_label, cumul=cumul)

        labels = self._get_labels(compare_label)

        if labels is None:
            labels = self.config_dict[compare_label]

        for i, f in enumerate(self.f_list):
            # print(len(self.config_dict[compare_label]))
            # if i < len(self.config_dict[compare_label]):
            sns.distplot(f,
                         hist=False,
                         kde_kws={
                             'shade': True,
                             'lw': 2,
                             'zorder': i
                         },
                         color=self.color_list[i],
                         label=labels[i])

        if self.sim:
            f_true = Results(self.sim_output_file).get_truths(['f'])['f']
            ax.axvline(f_true, 0, 10, color='k', zorder=3, lw=2., alpha=0.7)

        self._annotate(ax, title=title, extend=extend)

        fig.savefig(self.savefile_templ, bbox_inches='tight', dpi=300)

    def _annotate(self, ax, title=None, extend=False):
        '''Plot the annotations'''
        # annotations
        ax.set_xlim(0, 1)
        if title is not None:
            ax.set_title(title, fontsize=24)
        ax.set_xlabel('$f$')
        ax.set_ylabel('$P(f | \hat{E}, \hat{\omega})$')

        legend_loc = (1.8, 1) if extend else (1.3, 1)
        ax.legend(fontsize=22, bbox_to_anchor=legend_loc)