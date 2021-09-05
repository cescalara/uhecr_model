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
output_path = os.path.join(path_to_this_file, "..", "output")


class OutputFigures():
    def __init__(self,
                 fig_args,
                 sim_model=None,
                 sim_output_file=None,
                 tight_B=False,
                 skymap_label_fmt="default"):
        '''Collection of figures obtained from output of fitting.'''
        # parameters set from argparse
        self.source_type = fig_args["source"]
        self.detector_type = fig_args["detector"]
        self.model_type = fig_args["model"]
        self.data_type = fig_args["dtype"]
        self.ptype = fig_args["ptype"]
        self.seed = fig_args["seed"]
        self.sim_model_type = sim_model
        self.sim_output_file = sim_output_file

        # output file from simulation / fitting
        if tight_B:
            self.output_file = os.path.join(
                output_path,
                "tmp_{0}_fit_{5}_{1}_{2}_{3}_{4}_tightB.h5".format(
                    self.model_type, self.source_type, self.detector_type,
                    self.seed, self.ptype, self.data_type))
        else:
            self.output_file = os.path.join(
                output_path, "tmp_{0}_fit_{5}_{1}_{2}_{3}_{4}.h5".format(
                    self.model_type, self.source_type, self.detector_type,
                    self.seed, self.ptype, self.data_type))

        if sim_output_file is not None:
            self.output_file = os.path.join(
                output_path, "{0}_fit_{5}_{1}_{2}_{3}_{4}_{6}.h5".format(
                    self.model_type, self.source_type, self.detector_type,
                    self.seed, self.ptype, self.data_type,
                    self.sim_model_type))
            # if self.detector_type == "joint_gmf":
            #     self.output_file = os.path.join(
            #         output_path, "{0}_fit_{5}_{1}_{2}_{3}_{4}_{5}.h5".format(
            #             self.model_type, self.source_type, self.detector_type,
            #             self.seed, self.ptype, self.data_type,
            #             self.sim_model_type))

            # else:
            #     self.output_file = os.path.join(
            #         output_path, "{0}_fit_{5}_{1}_{2}_{3}_{4}.h5".format(
            #             self.model_type, self.source_type, self.detector_type,
            #             self.seed, self.ptype, self.data_type))

        # obtain detector properties / params from imports
        self.detector_properties, self.detector_params = self._get_detectorimports(
        )

        # initialize the data object and add relevant data into it
        self._initialize_data()

        # format style for skymap labels (from "mpl" or "TA")
        self.skymap_label_fmt = skymap_label_fmt

    def _get_detectorimports(self):
        '''Get variables imported by (detector_name).py'''
        if self.detector_type == "TA2015":
            from fancy.detector.TA2015 import detector_properties, detector_params
        elif self.detector_type == "auger2014":
            from fancy.detector.auger2014 import detector_properties, detector_params
        elif self.detector_type == "auger2010":
            from fancy.detector.auger2010 import detector_properties, detector_params
        else:
            raise Exception("Undefined detector type!")

        return detector_properties, detector_params

    def _initialize_data(self):
        '''Initialize data object based on sim / real data'''
        self.data = Data()

        if self.data_type == "sim":
            self.data.from_file(self.sim_output_file)

        elif self.data_type == "real":
            self.data.add_source(source_file, self.source_type)
            self.data.add_uhecr(uhecr_file, self.detector_type)
            self.data.add_detector(self.detector_properties)

    def src_uhecr_skymap(self, savefile=None, coord="G", exposure="map"):
        '''
        Plot skymap with sources + UHECR, obtained from Data object.
        Optionally add either exposure map or exposure limit.
        '''
        # labels
        src_label = self.data.source.label
        uhecr_label = self.data.uhecr.label
        # read in initializations
        omega_src = Direction(self.data.source.unit_vector)
        omega_arr = Direction(self.data.uhecr.unit_vector)
        energy_arr = self.data.uhecr.energy

        # convert source / arrival directions to equatorial / galactic
        if coord == "G":
            x_src, y_src = omega_src.glons, omega_src.glats
            x_arr, y_arr = omega_arr.glons, omega_arr.glats
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
        skymap.title("{0} + {1}".format(src_label, uhecr_label))

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

        skymap.save(savefile, bbox_inches='tight')

    def eval_association_probs(self):
        '''
        Evaluate the assocation probabilities between sources and UHECR
        from some output file
        '''

        # Log probability
        results = Results(self.output_file)
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

        # sort so that those with largest assossations appear first
        N_assos = {}
        for i in inds:
            N_assos[self.data.source.name[i].decode("UTF-8")] = len(
                np.argwhere([d == i for d in dominant]))

        N_assos_sorted = {
            k: v
            for k, v in sorted(N_assos.items(), key=lambda item: item[1])[::-1]
        }

        # print("Dominant sources: ", [self.data.source.name[i] for i in inds])

        return uhecr_p, N_assos_sorted

    def association_skymap(self, savefile, coord="G"):
        '''Plot association skymap between sources and UHECRs'''

        uhecr_p, N_assos_sorted = self.eval_association_probs()

        # labels
        # src_label = self.data.source.label

        src_label = "Swift-BAT" if self.data.source.label == "swift_BAT_213" else self.data.source.label.split(
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
        skymap.save(savefile)

        return N_assos_sorted

    def corner(self, savefile):
        '''Plot corner plot'''

        if self.dtype == "sim":
            self.corner_sim(savefile)
        elif self.dtype == "real":
            self.corner_data(savefile)

    def corner_data(self, savefile):
        '''Plot corner plot of data'''
        # Get chains from joint fit
        results_fit = Results(self.output_file)

        # get keys and corresponding labels
        labels = {}
        if self.model_type.find("join") != -1:
            keys = ['alpha', 'B', 'f']

            labels['B'] = r'$B$ / $\mathrm{nG}$'
            labels['alpha'] = r'$\alpha$'
            labels['f'] = r'$f$'

        elif self.model_type == "arrival_direction":
            keys = ['kappa', 'L', 'f']

            labels['L'] = r'$L$'
            labels['kappa'] = r'$\kappa$'
            labels['f'] = r'$f$'

        chain = results_fit.get_chain(keys)

        # Make nicely labelled dict
        chain_for_df = {}
        for key in keys:
            chain_for_df[labels[key]] = chain[key]

        # Make ordered dataframe
        df = DataFrame(data=chain_for_df)
        df = df[[labels[keys[0]], labels[keys[1]], labels[keys[2]]]]

        corner = Corner(df, color=midblue, contour_color=midblue_contour)
        corner.save(savefile)

    def corner_sim(self, savefile, cumul=False, output_files=None):

        if cumul:
            self.corner_sim_cumul(output_files, savefile)
        else:
            results_fit = Results(self.output_file)

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

            if self.model_type == "arrival_direction":
                self.corner_sim_arrival(results_fit, truth, savefile)
            else:
                self.corner_sim_joint(results_fit, truth, savefile)

    def corner_sim_arrival(self, results_fit, truth, savefile):
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
                        contour_color=midblue_contour)

        corner.save(savefile)

    def corner_sim_joint(self, results_fit, truth, savefile):
        '''Plot corner plot resulting from simulations'''

        keys = ['F0', 'L', 'alpha', 'B', 'f']
        chain = results_fit.get_chain(keys)

        # Convert form Stan units to plot units
        chain['F0'] = chain['F0'] / 1.0e3  # km^-2 yr^-1
        chain['L'] = chain['L'] * 10  # 10^-38 yr^-1

        labels = {}
        labels['L'] = r'$L$ / $10^{38}$ $\mathrm{yr}^{-1}$'
        labels['F0'] = r'$F_0$ / $\mathrm{km}^{-2} \ \mathrm{yr}^{-1}$'
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
        df = df[[
            labels['F0'], labels['L'], labels['alpha'], labels['B'],
            labels['f']
        ]]

        corner = Corner(df, truths, color=purple, contour_color=purple_contour)
        corner.save(savefile)

    def corner_sim_cumul(self, output_files, savefile):
        '''Cumulative distribution of corner plots from simulation'''

        keys = ['F0', 'L', 'alpha', 'B', 'f']
        chain_avgs = {key: 0 for key in keys}
        chain_list = []
        Nseeds = len(output_files)

        # get chains for each output file
        # ignores self.outpuf_file
        for output_file in output_files:
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
        chain_avgs['F0'] = chain_avgs['F0'] / 1.0e3  # km^-2 yr^-1
        chain_avgs['L'] = chain_avgs['L'] * 10  # 10^-38 yr^-1

        # Get truths from simulation
        results_sim = Results(self.sim_output_file)

        truth_keys = ['F0', 'L', 'alpha', 'B', 'f']
        truth = results_sim.get_truths(truth_keys)
        info_keys = ['Eth', 'Eth_sim']
        info = results_sim.get_truths(info_keys)

        # Correct for different Eth in sim and fit
        # Also scale to plot units
        flux_scale = (info['Eth'] / info['Eth_sim'])**(1 - truth['alpha'])
        truth['F0'] = truth['F0'] * flux_scale  # km^-2 yr^-1
        truth['L'] = truth['L'][0] * flux_scale / 1.0e39 * 10  # 10^-38 yr^-1

        labels = {}
        labels['L'] = r'$L$ / $10^{38}$ $\mathrm{yr}^{-1}$'
        labels['F0'] = r'$F_0$ / $\mathrm{km}^{-2} \ \mathrm{yr}^{-1}$'
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
        df = df[[
            labels['F0'], labels['L'], labels['alpha'], labels['B'],
            labels['f']
        ]]

        corner = Corner(df, truths, color=purple, contour_color=purple_contour)
        corner.save(savefile)


class SourceUHECRDist():
    def __init__(self, figsize=(6, 4), sim_output_file=None):
        '''Class that organizes plotting of source-UHECR association fraction distribution'''
        self.fig, self.ax = plt.subplots(figsize=figsize)

        self.color_list = [grey, lightpurple, lightblue]

        self.sim_output_file = sim_output_file

    def get_fs(self, fname_list, cumul=False):
        '''Get source-UHECR association fraction from output files'''

        self.f_list = []

        if cumul:  # evaluate cumulative distribution
            # then each element in fname_list must be a list of fnames with different models
            # each of this list will contain fnames with different seeds
            for fname_models_list in fname_list:
                f_list = []
                for output_file in fname_models_list:
                    f_i = Results(output_file).get_chain(['f'])['f']
                    f_list.append(f_i)

                f_avg = np.mean(np.array(f_list), axis=0)
                self.f_list.append(f_avg)

        else:
            for output_file in fname_list:
                try:
                    f_i = Results(output_file).get_chain(['f'])['f']
                except KeyError:
                    f_i = np.zeros(100)
                self.f_list.append(f_i)

        return self.f_list

    def plotdist(self, labels, extend=False, title=None):
        '''Plot the distribution'''
        for i, label in enumerate(labels):
            # for arrival_direction -> arrival direction
            label = label.replace("_", " ") if "_" in label else label

            sns.distplot(self.f_list[i],
                         hist=False,
                         kde_kws={
                             'shade': True,
                             'lw': 2,
                             'zorder': i
                         },
                         color=self.color_list[i],
                         label=label)

        if self.sim_output_file is not None:
            f_true = Results(self.sim_output_file).get_truths(['f'])['f']
            self.ax.axvline(f_true,
                            0,
                            10,
                            color='k',
                            zorder=3,
                            lw=2.,
                            alpha=0.7)

        self._annotate(title=title, extend=extend)

    def _annotate(self, title=None, extend=False):
        '''Plot the annotations'''
        # annotations
        self.ax.set_xlim(0, 1)
        if title:
            self.ax.set_title(title, fontsize=24)
        self.ax.set_xlabel('$f$')
        self.ax.set_ylabel('$P(f | \hat{E}, \hat{\omega})$')

        legend_loc = (1.6, 1) if extend else (1.3, 1)
        self.ax.legend(fontsize=22, bbox_to_anchor=legend_loc)

    def save(self, savefile):
        '''Save the distribution'''
        # save
        self.fig.savefig(savefile, bbox_inches='tight')