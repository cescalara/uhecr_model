'''Collection of functions that generate skymaps related to simulation of arrival directions with consideration to GMF deflections.'''

import healpy
import os
import h5py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
from scipy.interpolate import griddata

from fancy.plotting.allskymap_cartopy import AllSkyMapCartopy as AllSkyMap
from fancy.interfaces.stan import uv_to_coord


class GMFSimFigures():
    def __init__(self, sim_output_file):
        '''Container for plots resulting from GMF simulations'''

        self.sim_output_file = sim_output_file
        # read all plotting variables from simulation output
        self._read_from_file()

    def _read_from_file(self):
        '''Read plotting variables from simulation output file.'''
        with h5py.File(self.sim_output_file, "r") as f:

            self.src_indices = f["plotvars/src_indices"][()]
            self.bg_indices = f["plotvars/bg_indices"][()]

            self.coord_gb = uv_to_coord(f["plotvars/omega_gb"][()])
            self.coord_earth = uv_to_coord(f["plotvars/omega_earth"][()])
            self.coord_det_exp_limited = uv_to_coord(
                f["plotvars/omega_det_exp_limited"][()])
            self.coord_det = uv_to_coord(f["plotvars/omega_det"][()])

            self.energies_earth = f["plotvars/energies_earth"][()]
            self.energies_det_exp_limited = f[
                "plotvars/energies_det_exp_limited"][()]
            self.energies_det = f["plotvars/energies_det"][()]

            self.map_unlensed = f["plotvars/map_unlensed"][()]
            self.map_lensed = f["plotvars/map_lensed"][()]
            self.energies_nlarge = f["plotvars/energies_nlarge"][()]
            self.coord_nlarge = uv_to_coord(f["plotvars/omega_nlarge"][()])

            self.coord_src = uv_to_coord(f["source/unit_vector"][()])

            detector_label = f["detector/label"][()].decode("UTF-8")

        # get detector params for plotting exposure map
        self._get_detector_params(detector_label)

    def _get_detector_params(self, detector_type):
        '''set detector and detector params'''
        if detector_type == "TA2015":
            from fancy.detector.TA2015 import detector_params
        elif detector_type == "auger2014":
            from fancy.detector.auger2014 import detector_params
        elif detector_type == "auger2010":
            from fancy.detector.auger2010 import detector_params
        else:
            raise Exception("Undefined detector type!")

        self.detector_params = detector_params

    def uhecr_gb_skymap(self, savefile, trunc=True):
        '''
        Plot skymap containing UHECRs randomly sampled from vMF / spherical
        distribution using stan code, up until the galactic boundary.
        '''
        lons_gb, lats_gb = 180 - self.coord_gb.galactic.l.deg, self.coord_gb.galactic.b.deg
        lons_src, lats_src = 180 - self.coord_src.galactic.l.deg, self.coord_src.galactic.b.deg

        cutoff_idx = 50 if trunc else 100  # plot only half of the random points (Nrand=100)

        skymap = AllSkyMap(lon_0=180)
        skymap.set_gridlines(color="k", linestyle=":", lw=10, zorder=1)
        skymap.scatter(lons_gb[:cutoff_idx, self.src_indices],
                       lats_gb[:cutoff_idx, self.src_indices],
                       s=4.5,
                       color="r",
                       alpha=0.25,
                       lw=0.,
                       label="UHECR from Sources")
        skymap.scatter(lons_gb[:cutoff_idx, self.bg_indices],
                       lats_gb[:cutoff_idx, self.bg_indices],
                       s=4.5,
                       color="b",
                       alpha=0.25,
                       lw=0.,
                       label="UHECR from Background")
        skymap.scatter(lons_src,
                       lats_src,
                       s=60.0,
                       color="g",
                       alpha=0.9,
                       marker="*",
                       label="Sources")

        handles = [
            mlines.Line2D([], [],
                          color='r',
                          marker='o',
                          lw=0,
                          markersize=4,
                          alpha=0.3,
                          label="UHECR from Sources"),
            mlines.Line2D([], [],
                          color='b',
                          marker='o',
                          lw=0,
                          markersize=4,
                          alpha=0.3,
                          label="UHECR from Background"),
            mlines.Line2D([], [],
                          color='g',
                          marker='*',
                          lw=0,
                          markersize=10,
                          label="Sources")
        ]

        skymap.legend(handles=handles, bbox_to_anchor=(0.9, 0.2))
        skymap.title(
            "Skymap of vMF distributed UHECRs from Sources and Background")

        skymap.save(savefile)

    def _interpolate_crmap(self, crmap, NPIX=49152):
        '''
        Perform 2-D interpolation on healpix map to allow plotting with matplotlib.
        
        NPIX is set via ParticleMapContainer in CRPropa.
        '''
        # get coordinates from pixels
        ipixs = np.arange(0, NPIX, 1, dtype=int)
        th, ph = healpy.pix2ang(healpy.npix2nside(NPIX), ipixs)
        lons_crmap, lats_crmap = -np.rad2deg(np.pi - ph), np.rad2deg(np.pi /
                                                                     2. - th)

        max_val = np.max(crmap)  # to normalize

        grid_lons, grid_lats = np.mgrid[-180:180:300j, -90:90:300j]

        # interpolation
        grid_crmap = griddata((lons_crmap, lats_crmap),
                              crmap / max_val, (grid_lons, grid_lats),
                              method="cubic",
                              fill_value=0)

        return grid_lons, grid_lats, grid_crmap

    def plot_crmap(self, lons, lats, crmap, title, savefile):
        '''Plot UHECR probability map'''

        skymap = AllSkyMap(lon_0=180)
        skymap.set_gridlines(color="k", linestyle=":", lw=10, zorder=1)

        im = skymap.contourf(lons,
                             lats,
                             crmap,
                             cmap=cm.afmhot_r,
                             levels=np.linspace(0, 1, 50),
                             alpha=1)
        skymap.ax.set_facecolor('w')
        skymap.title(title)
        cbar = skymap.fig.colorbar(im,
                                   pad=0.1,
                                   cmap=cm.afmhot_r,
                                   orientation='horizontal',
                                   alpha=1,
                                   shrink=0.65)

        cbar.set_ticks(np.arange(0, 1.1, 0.1))
        cbar.set_ticklabels(
            ["{0:.1f}".format(i) for i in np.arange(0, 1.1, 0.1)])

        skymap.save(savefile)

    def unlensed_skymap(self, savefile, healpix=False):
        '''
        Plot unlensed skymap of UHECR at galactic boundary using ax.contourf. 
        Normalized from [0, 1] to represent probability map.
        '''
        title = "UHECR Map at Galactic Boundary"
        if healpix:  # view using healpix map
            healpy.mollview(map=self.map_unlensed, title=title, cmap=cm.hot)
        else:  # view using matplotlib contourf + griddata interpolation
            grid_lons, grid_lats, grid_crmap_unlensed = self._interpolate_crmap(
                self.map_unlensed)
            self.plot_crmap(grid_lons, grid_lats, grid_crmap_unlensed, title,
                            savefile)

    def lensed_skymap(self, savefile, healpix=False):
        '''
        Plot lensed skymap of UHECR at earth using ax.contourf. 
        Normalized from [0, 1] to represent probability map.
        '''
        title = "UHECR Map at Earth"
        if healpix:  # view using healpix map
            healpy.mollview(map=self.map_lensed, title=title, cmap=cm.hot)
        else:  # view using matplotlib contourf + griddata interpolation
            grid_lons, grid_lats, grid_crmap_lensed = self._interpolate_crmap(
                self.map_lensed)
            self.plot_crmap(grid_lons, grid_lats, grid_crmap_lensed, title,
                            savefile)

    def sim_uhecr_skymap(self, savefile, exposure=True, trunc=True):
        '''
        Plot resulting simulated UHECR from GMF simulation onto skymap. Exposure = True will
        plot the exposure map onto the plot.
        '''
        lons_det, lats_det = 180 - self.coord_det.galactic.l.deg, self.coord_det.galactic.b.deg
        lons_earth, lats_earth = 180 - self.coord_nlarge.galactic.l.deg, self.coord_nlarge.galactic.b.deg
        lons_src, lats_src = 180 - self.coord_src.galactic.l.deg, self.coord_src.galactic.b.deg

        skymap = AllSkyMap(lon_0=180)
        skymap.set_gridlines(color="k", linestyle=":", lw=10, zorder=1)

        cutoff_idx = 13400 if trunc else 26800  # plot only half (Nuhecr*Nrand = 26800)

        skymap.scatter(lons_det,
                       lats_det,
                       color="r",
                       alpha=0.5,
                       label="UHECR_truncated",
                       marker="o",
                       zorder=5,
                       lw=0)
        skymap.scatter(lons_earth[:cutoff_idx],
                       lats_earth[:cutoff_idx],
                       s=2.0,
                       color="b",
                       alpha=0.075,
                       label="UHECR_Earth",
                       marker="o",
                       zorder=4,
                       lw=0)
        skymap.scatter(lons_src,
                       lats_src,
                       s=60.0,
                       color="g",
                       alpha=1.0,
                       marker="*",
                       label="Sources",
                       zorder=10)

        if exposure:
            skymap.exposure_map(self.detector_params, zorder=1)

        handles = [
            mlines.Line2D([], [],
                          color='r',
                          marker='o',
                          lw=0,
                          markersize=8,
                          alpha=0.3,
                          label="Simulated UHECR"),
            mlines.Line2D([], [],
                          color='b',
                          marker='o',
                          lw=0,
                          markersize=2,
                          alpha=0.3,
                          label="Sampled UHECR"),
            mlines.Line2D([], [],
                          color='g',
                          marker='*',
                          lw=0,
                          markersize=10,
                          label="Sources")
        ]

        skymap.legend(handles=handles, bbox_to_anchor=(0.9, 0.2))
        skymap.title("Skymap of Simulated UHECR for TA Hotspot")
        skymap.save(savefile)