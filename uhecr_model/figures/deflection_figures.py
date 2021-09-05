'''Collection of functions that generate deflection skymaps.'''

import matplotlib.pyplot as plt
import numpy as np
import os
import h5py
import matplotlib.lines as mlines
from astropy.coordinates import SkyCoord
from astropy import units as u

from fancy import Uhecr
from fancy.interfaces.stan import uv_to_coord
from fancy.plotting.allskymap_cartopy import AllSkyMapCartopy as AllSkyMap
from fancy.interfaces.stan import coord_to_uv
from scipy.optimize import root

# use minimalist style
plt.style.use("minimalist")

# paths to important files
path_to_this_file = os.path.abspath(os.path.dirname(__file__))
stan_path = os.path.join(path_to_this_file, "..", "stan")
source_file = os.path.join(path_to_this_file, "..", "data", "sourcedata.h5")
uhecr_file = os.path.join(path_to_this_file, "..", "data", "UHECRdata.h5")
table_path = os.path.join(path_to_this_file, "..", "tables")
output_path = os.path.join(path_to_this_file, "..", "output")
plot_path = os.path.join(path_to_this_file, "..", "..", "..", "plots")


class DeflectionFigures():
    def __init__(self, detector="TA2015", header="tmp", gmf_model="JF12"):
        '''Figures that demonstrate GMF deflections and evaluation of kappa_gmf'''
        self.ptypes_list = ["p", "He", "N", "O", "Si", "Fe"]
        self.gmf_model = gmf_model
        # read kappa_gmf values for each particle type
        self.read_from_file(detector)

        self.header = header

    def read_from_file(self, detector):
        '''Read information needed for plotting from the UHECR file'''

        self.kappa_gmf_list = []
        self.coord_gal_list = []
        self.coord_rand_list = []

        self.coord_true_list = []
        # note:
        # omega_true == uhecr lons + lats (N_uhecr, 2)
        # omega_rand == vMF sampled lons + lats (N_uhecr, Nrand, 2)
        # omega_gal == deflected lons + lats (N_uhecr, Nrand, 2)

        with h5py.File(uhecr_file, "r") as f:

            detector_dict = f[detector]

            kappa_gmf_dict = detector_dict["kappa_gmf"]
            gmf_model_dict = kappa_gmf_dict[self.gmf_model]

            for ptype in self.ptypes_list:
                kappa_gmf_dict_ptype = gmf_model_dict[ptype]
                self.coord_true_list.append(
                    np.rad2deg(kappa_gmf_dict_ptype["omega_true"][()]))
                self.coord_rand_list.append(
                    np.rad2deg(kappa_gmf_dict_ptype["omega_rand"][()]))
                self.coord_gal_list.append(
                    np.rad2deg(kappa_gmf_dict_ptype["omega_gal"][()]))
                self.kappa_gmf_list.append(
                    kappa_gmf_dict_ptype["kappa_gmf"][()])

        self.N_uhecr = int(self.coord_rand_list[0].shape[0])
        self.Nrand = int(self.coord_rand_list[0].shape[1])

        # self.min_energy_idx = int(np.argmin(detector_dict["energy"][()]))
        self.min_energy_idx = 3

    def _create_savefile(self, ptype, sel_uhecr=False, nodefl=False):
        '''Create savefile based on given conditions'''

        savefile_templ = "{0}_{1}_{2}_{{0}}.png".format(
            self.header, ptype, self.gmf_model)
        if sel_uhecr is not None:
            savefile_templ = "{0}_{1}_{2}_{{0}}.png".format(
                self.header, ptype, sel_uhecr)
        if nodefl:
            savefile_templ = "{0}_{1}_nodefl_{{0}}.png".format(
                self.header, ptype)

        # create header directories if they dont exist
        header_path = os.path.join(plot_path, self.header)
        if not os.path.exists(header_path):
            os.mkdir(header_path)

        if not os.path.exists(os.path.join(header_path, "defl")):
            os.mkdir(os.path.join(header_path, "defl"))

        savefile_templ = os.path.join(header_path, "defl", savefile_templ)

        return savefile_templ

    def plot_smearing(self, ptype="p", sel_uhecr_idx=None, nodefl=False):
        '''Plot the smeared out deflections onto skymap'''

        savefile_templ = self._create_savefile(ptype, sel_uhecr_idx, nodefl)

        skymap = AllSkyMap(lon_0=180)
        skymap.set_gridlines(label_fmt="default")
        # get corresponding index where the particle desired is
        ptype_idx = np.argwhere([p == ptype for p in self.ptypes_list])[0][0]
        # title = "Deflection skymap with TA hotspot data - {0}".format(ptype)

        # plot all of them as a scatter plot
        skymap.scatter(
            self.coord_true_list[ptype_idx][:, 0],  # lons_true
            self.coord_true_list[ptype_idx][:, 1],  # lats_true
            color="k",
            alpha=1,
            marker="+",
            s=20.0)
        skymap.scatter(
            self.coord_rand_list[ptype_idx][:, :, 0],  # lons_rand
            self.coord_rand_list[ptype_idx][:, :, 1],  # lats_rand
            color="b",
            alpha=0.1,
            s=5.0,
            lw=0)

        if not nodefl:
            skymap.scatter(
                self.coord_gal_list[ptype_idx][:, :, 0],  # lons_defl
                self.coord_gal_list[ptype_idx][:, :, 1],  # lats_defl
                color="r",
                alpha=0.1,
                s=5.0,
                lw=0)

        handles = [
            mlines.Line2D([], [],
                          color='k',
                          marker='+',
                          lw=0,
                          markersize=6,
                          alpha=1,
                          label="Original UHECR"),
            mlines.Line2D([], [],
                          color='b',
                          marker='o',
                          lw=0,
                          markersize=3,
                          alpha=0.3,
                          label="Randomized UHECR"),
            mlines.Line2D([], [],
                          color='r',
                          marker='o',
                          lw=0,
                          markersize=3,
                          alpha=0.3,
                          label="Deflected UHECR")
        ]

        if nodefl:
            # title = "vMF-sampled skymap using $\sigma_\omega$ with TA hotspot data"
            handles = handles[:-1]

        legend1 = skymap.legend(handles=handles,
                                bbox_to_anchor=(1.25, 0.27),
                                fontsize=14)
        skymap.ax.add_artist(legend1)

        # if we want to highlight a single particle
        if sel_uhecr_idx is not None:

            # select that with minimum energy if requested
            if sel_uhecr_idx == "Emin":
                sel_uhecr_idx = self.min_energy_idx  # UHECR with lowerst energy
            else:
                sel_uhecr_idx = sel_uhecr_idx
            sel_uhecr_lonlat = 180. - self.coord_true_list[ptype_idx][
                sel_uhecr_idx,
                0], self.coord_true_list[ptype_idx][sel_uhecr_idx, 1]

            # rand -> defl for a single UHECR
            skymap.scatter(self.coord_rand_list[ptype_idx][sel_uhecr_idx, :,
                                                           0],
                           self.coord_rand_list[ptype_idx][sel_uhecr_idx, :,
                                                           1],
                           color="g",
                           alpha=0.1,
                           s=5.0,
                           lw=0)
            skymap.scatter(self.coord_gal_list[ptype_idx][sel_uhecr_idx, :, 0],
                           self.coord_gal_list[ptype_idx][sel_uhecr_idx, :, 1],
                           color="g",
                           alpha=0.8,
                           s=5.0,
                           lw=0)

            # plot lines between each random direction and deflected one
            for i in range(self.Nrand):
                skymap.geodesic(self.coord_rand_list[ptype_idx][sel_uhecr_idx,
                                                                i, 0],
                                self.coord_rand_list[ptype_idx][sel_uhecr_idx,
                                                                i, 1],
                                self.coord_gal_list[ptype_idx][sel_uhecr_idx,
                                                               i, 0],
                                self.coord_gal_list[ptype_idx][sel_uhecr_idx,
                                                               i, 1],
                                alpha=0.05,
                                color="k")

            # sel_uhecr_handles = [
            #     mlines.Line2D(
            #         [], [],
            #         color='g',
            #         marker='o',
            #         lw=0,
            #         markersize=4,
            #         alpha=0.3,
            #         label=
            #         "UHECR Coordinate (Gal):\n ({0:.1f}$^\circ$, {1:.1f}$^\circ$)"
            #         .format(*sel_uhecr_lonlat))
            # ]

            # legend2 = skymap.legend(handles=sel_uhecr_handles,
            #                         bbox_to_anchor=(1.25, 1.05))

            # skymap.ax.add_artist(legend2)

        # skymap.title(title)
        skymap.save(savefile_templ.format("smear_skymap"))

    def plot_circles(self, ptype="p", sel_uhecr_idx=None, ang_err=1.7):
        '''
        Plot the circles that represent the effective smearing applied to
        our model onto the skymap.
        '''

        savefile_templ = self._create_savefile(ptype,
                                               sel_uhecr_idx,
                                               nodefl=False)

        ptype_idx = np.argwhere([p == ptype for p in self.ptypes_list])[0][0]
        # print(kappa_gmf_list[ptype_idx])

        # evaluate the deflection angle that corresponds to the radius of
        # deflection circle
        cos_theta_arr = np.zeros(len(self.kappa_gmf_list[ptype_idx]))
        for i, kappa_gmf in enumerate(self.kappa_gmf_list[ptype_idx]):
            sol = root(fischer_int_eq_P, x0=1, args=(kappa_gmf, 0.683))
            cos_theta = sol.x[0]
            cos_theta_arr[i] = cos_theta

        thetas = np.rad2deg(np.arccos(cos_theta_arr))

        skymap = AllSkyMap(lon_0=180)
        skymap.set_gridlines(label_fmt="default")

        skymap.scatter(self.coord_true_list[ptype_idx][:, 0],
                       self.coord_true_list[ptype_idx][:, 1],
                       color="k",
                       alpha=1,
                       marker="+",
                       s=10.0)

        for i, (lon, lat) in enumerate(self.coord_true_list[ptype_idx][:]):
            skymap.tissot(lon, lat, ang_err, color="b", alpha=0.2)
            skymap.tissot(lon, lat, thetas[i], color="r", alpha=0.05)

        handles = [
            mlines.Line2D([], [],
                          color='k',
                          marker='+',
                          lw=0,
                          markersize=8,
                          alpha=1,
                          label="Original UHECR"),
            mlines.Line2D([], [],
                          color='b',
                          marker='o',
                          lw=0,
                          markersize=8,
                          alpha=0.3,
                          label=r"$\sigma_\omega$"),
            mlines.Line2D([], [],
                          color='r',
                          marker='o',
                          lw=0,
                          markersize=8,
                          alpha=0.3,
                          label=r"$\sigma_{{\omega + GMF}}$")
        ]

        # if we want to highlight a single particle
        if sel_uhecr_idx is not None:

            # select that with minimum energy if requested
            if sel_uhecr_idx == "Emin":
                sel_uhecr_idx = self.min_energy_idx  # UHECR with lowerst energy
            else:
                sel_uhecr_idx = sel_uhecr_idx
            sel_uhecr_lonlat = 180. - self.coord_true_list[ptype_idx][
                sel_uhecr_idx,
                0], self.coord_true_list[ptype_idx][sel_uhecr_idx, 1]

            skymap.tissot(self.coord_true_list[ptype_idx][sel_uhecr_idx, 0],
                          self.coord_true_list[ptype_idx][sel_uhecr_idx, 1],
                          ang_err,
                          color="b",
                          alpha=0.2)
            skymap.tissot(self.coord_true_list[ptype_idx][sel_uhecr_idx, 0],
                          self.coord_true_list[ptype_idx][sel_uhecr_idx, 1],
                          thetas[sel_uhecr_idx],
                          color="g",
                          alpha=0.3)

        skymap.legend(handles=handles, bbox_to_anchor=(1.3, 0.27), fontsize=14)
        # skymap.title("Deflection skymap with data - {0}".format(ptype))
        skymap.save(savefile_templ.format("circle_skymap"))

    def histogram(self, ptypes=["p", "N", "Fe"]):
        '''Create histogram of kappa_gmf'''

        savefile_templ = self._create_savefile(ptype="p_N_Fe")

        fig, ax = plt.subplots(figsize=(6, 4))

        for ptype in ptypes:
            ptype_idx = np.argwhere([p == ptype
                                     for p in self.ptypes_list])[0][0]

            ax.hist(self.kappa_gmf_list[ptype_idx])

        fig.savefig(savefile_templ.format("histogram"))
        pass


def fischer_int(kappa, cos_thetaP):
    '''Integral of vMF function over all angles'''
    return (1. - np.exp(-kappa *
                        (1 - cos_thetaP))) / (1. - np.exp(-2. * kappa))


def fischer_int_eq_P(cos_thetaP, kappa, P):
    '''Equation to find roots for'''
    return fischer_int(kappa, cos_thetaP) - P