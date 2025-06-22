# This file is part of cp_verify.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This software is dual licensed under the GNU General Public License and also
# under a 3-clause BSD license. Recipients may choose which of these licenses
# to use; please see the files gpl-3.0.txt and/or bsd_license.txt,
# respectively.  If you choose the GPL option then the following text applies
# (but note that there is still no warranty even if you opt for BSD instead):
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import os
import argparse
import logging

import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from matplotlib.pyplot import cm
from matplotlib.colors import SymLogNorm

from lsst.daf.butler import Butler
from lsst.daf.butler import EmptyQueryResultError, DatasetNotFoundError

from .utils import plot_fp_statistic


class DetailedPlotter():
    """A class to generate calibration verification reports.

    Parameters
    ----------
    repo : `str`
        The location of the butler repository to retrieve results
        from.
    instrument : `str`
        The instrument associated with this data.
    output_path : `str`
        The location the report will be written to.
    collections : `list` [`str`]
        A list of collections to search.
    **kwargs :
        Other keyword parameters.  Currently parsed values:

        - ``do_copy``: Should files be copied from butler (`bool`).
        - ``do_overwrite``: Should pre-existing files be overwritten (`bool`).
    """

    def __init__(self, repo, output, collections, instrument, detectorIds,
                 plotSummaryStats, doOverwrite, **kwargs):
        super().__init__()
        # Set source and destination information.
        self.repo = repo
        self.output = os.path.join(output, "images")
        self.collections = collections
        self.instrument = instrument
        self.doOverwrite = doOverwrite
        self.doPlotSummaryStats = plotSummaryStats

        # Logging
        self.log = logging.getLogger(__name__) if "log" not in kwargs else kwargs["log"]

        # Instantiate a butler for our repository.
        self.butler = Butler(repo, collections=self.collections)
        self.registry = self.butler.registry

        # Get the camera associated with these calibrations
        self.camera = self.butler.get(
            "camera",
            instrument=self.instrument,
            collections=f"{self.instrument}/calib",
        )

        # Get the detector ids
        self.detectorIds = range(len(self.camera)) if len(detectorIds) == 0 else sorted(detectorIds)
        self.detectorIdsString = ','.join([str(_) for _ in self.detectorIds])

        # Start with no references to datasets and eventually
        # add the ones that exist.
        self.refs = None

        # Index the images so that they are saved in the file
        # system in a particular order. This is helpful for
        # navigation and formatting on the web interface
        self.plotIdx = 0

        # Common plotting parameters
        self.color = cm.tab20(np.linspace(0, 1, 16))

        # Collect all the data
        self.typesDict = self._makeEmptyAmpStatData()
        self.ptcTurnoffsDict = self._makeEmptyAmpStatData()
        self.linearityTurnoffsDict = self._makeEmptyAmpStatData()
        self.globalCtisDict = self._makeEmptyAmpStatData()
        self.serialCtiTurnoffsDict = self._makeEmptyAmpStatData()
        self.n00sDict = self._makeEmptyAmpStatData()
        self.empiricalN00sDict = self._makeEmptyAmpStatData()
        self.empiricalN00sResidualsDict = self._makeEmptyAmpStatData()
        self.n01sDict = self._makeEmptyAmpStatData()
        self.n10sDict = self._makeEmptyAmpStatData()
        self.gainsDict = self._makeEmptyAmpStatData()
        self.gainsUnadjustedDict = self._makeEmptyAmpStatData()
        self.gainsUnadjustedResidualsDict = self._makeEmptyAmpStatData()
        self.a00sDict = self._makeEmptyAmpStatData()

    def _getRaftBayFromRef(self, ref):
        detector = self.camera[ref.dataId['detector']]
        detName = detector.getName()
        raft = detName[0:3]
        bay = detName[4:]
        return raft, bay

    def _getRaftBayPathFromRef(self, ref):
        raft, bay = self._getRaftBayFromRef(ref)
        return f"{raft}/{bay}"

    def _getAllDictValues(self, d):
        def _nested(_):
            for v in _.values():
                if isinstance(v, dict):
                    yield from _nested(v)
                else:
                    yield v
        return np.asarray(list(_nested(d)))

    def _makeEmptyAmpStatData(self):
        empty = dict()
        for detector in self.camera:
            detName = detector.getName()
            empty[detName] = dict()
            for amp in detector:
                empty[detName][amp.getName()] = 0.0

        return empty

    def initializeOutputDirectories(self):
        self.log.info(f"Initializing output directories in {self.output}")
        if os.path.isdir(self.output) and len(os.listdir(self.output)) > 0:
            msg = "The output directory already exists and is not empty."
            if self.doOverwrite:
                self.log.warning(msg)
            else:
                raise RuntimeError(msg)

        # Initialize output directories
        for detector in self.camera:
            detName = detector.getName()
            raft = detName[0:3]
            bay = detName[4:]
            os.makedirs(f"{self.output}/{raft}", exist_ok=True)
            os.makedirs(f"{self.output}/{raft}/{bay}", exist_ok=True)

    def plotLinearizerStatistics(self):
        self.plotIdx += 1

        # Get datasets, if they exist
        try:
            refs = self.butler.query_datasets(
                "linearizer",
                instrument=self.instrument,
                where=f"detector in ({self.detectorIdsString})",
                collections=self.collections,
            )
        except EmptyQueryResultError:
            return

        self.refs = refs

        for ref in tqdm(self.refs, total=len(self.refs), desc="LINEARIZER"):
            lin = self.butler.get("linearizer", dataId=ref.dataId, collections=self.collections)
            ptc = self.butler.get("linearizerPtc", dataId=ref.dataId, collections=self.collections)

            fig, axes = plt.subplots(4, 4, figsize=(12, 12))
            fig.suptitle(self.camera[ref.dataId['detector']].getName() + f" ({ref.dataId['detector']})")
            axs = axes.ravel()
            detName = self.camera[ref.dataId['detector']].getName()
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                ampName = amp.getName()

                # Save dataset information
                self.typesDict[detName][ampName] = self.camera[ref.dataId['detector']].getPhysicalType()
                # Need to fix this. Why is it a zero dimensional array?
                linTurnoff = lin.linearityTurnoff[ampName].ravel()[0]
                self.linearityTurnoffsDict[detName][ampName] = linTurnoff

                # Plot
                ax = axs[i]
                ax.set_title(ampName)
                ax.scatter(ptc.photoCharges[ampName], ptc.rawMeans[ampName], s=1)
                ax.axhline(lin.linearityTurnoff[ampName], linestyle="--", color="k", alpha=0.5)

                ax.text(
                    2e-7, 0.25e5,
                    f"Linearity Turnoff:\n{round(linTurnoff, 1)} adu",
                )

                ax.set_xlim(-.5e-7, 6.5e-7)
                ax.set_ylim(0, 1.2e5)
                ax.ticklabel_format(style='sci', axis='both', scilimits=(0, 0), useMathText=True)
                ax.set_xlabel("Integrated PD\nCurrent")
                ax.set_ylabel("Raw Means (adu)")
            plt.tight_layout()
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}a-"
                        f"pd_current_and_raw_means_{self.camera[ref.dataId['detector']].getId():03}.png",
                        bbox_inches='tight')
            plt.close()

            fig, axes = plt.subplots(4, 4, figsize=(12, 12))
            fig.suptitle(self.camera[ref.dataId['detector']].getName() + f" ({ref.dataId['detector']})")
            axs = axes.ravel()
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                ampName = amp.getName()
                ax = axs[i]
                ax.set_title(ampName)
                ax.scatter(ptc.photoCharges[ampName], lin.fitResiduals[ampName]/ptc.rawMeans[ampName], s=1)
                ax.axhline(0, linestyle="--", color="k", alpha=0.5)
                ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)
                ax.set_ylim(-0.00075, 0.00075)
                if ref.dataId['detector'] in [60, 61, 62]:
                    ax.set_ylim(-0.003, 0.003)
                ax.set_xlabel("Integrated PD\nCurrent")
                ax.set_ylabel("fitResiduals / rawMeans")
            plt.tight_layout()
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}b-"
                        f"non-linearity_fit_residuals_{self.camera[ref.dataId['detector']].getId():03}.png",
                        bbox_inches='tight')
            plt.close()

    def plotPtcStatistics(self):
        self.plotIdx += 1

        # Get datasets, if they exist
        try:
            refs = self.butler.query_datasets(
                "ptc",
                instrument=self.instrument,
                where=f"detector in ({self.detectorIdsString})",
                collections=self.collections)
        except EmptyQueryResultError:
            return

        self.refs = refs

        for ref in tqdm(self.refs, total=len(self.refs), desc="PTC"):
            ptc = self.butler.get("ptc", dataId=ref.dataId, collections=self.collections)

            fig, axes = plt.subplots(4, 4, figsize=(12, 12))
            fig.suptitle(self.camera[ref.dataId['detector']].getName() + f" ({ref.dataId['detector']})")
            axs = axes.ravel()
            detName = self.camera[ref.dataId['detector']].getName()
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                ampName = amp.getName()

                # Save full array information
                self.typesDict[detName][ampName] = self.camera[ref.dataId['detector']].getPhysicalType()
                self.ptcTurnoffsDict[detName][ampName] = ptc.ptcTurnoff[ampName] * ptc.gain[ampName]
                self.n00sDict[detName][ampName] = ptc.noise[ampName]
                self.n01sDict[detName][ampName] = ptc.noiseMatrix[ampName][0][1]
                self.n10sDict[detName][ampName] = ptc.noiseMatrix[ampName][1][0]
                self.gainsDict[detName][ampName] = ptc.gain[ampName]
                self.gainsUnadjustedDict[detName][ampName] = ptc.gainUnadjusted[ampName]
                fracGainDiff = (ptc.gainUnadjusted[ampName] / ptc.gain[ampName]) - 1
                self.gainsUnadjustedResidualsDict[detName][ampName] = fracGainDiff
                self.a00sDict[detName][ampName] = ptc.aMatrix[ampName][0][0]
                empiricalNoise = np.nanmedian(ptc.noiseList[ampName])
                fracNoiseDiff = (empiricalNoise * ptc.gain[ampName]) / ptc.noise[ampName] - 1
                self.empiricalN00sResidualsDict[detName][ampName] = fracNoiseDiff
                self.empiricalN00sDict[detName][ampName] = fracNoiseDiff

                ax = axs[i]
                ax.set_title(ampName)
                ax.scatter(ptc.rawMeans[ampName], ptc.rawVars[ampName], s=1)
                ax.scatter(ptc.finalMeans[ampName], ptc.finalVars[ampName], s=1)

                if ampName in ptc.badAmps:
                    ax.set_facecolor('#fcd4d4')
                    ax.set_title(ampName + " (BAD)")

                ptcTurnoff = float(ptc.ptcTurnoff[ampName])
                gain = float(ptc.gain[ampName])
                a00 = float(ptc.aMatrix[ampName][0][0])
                n00 = float(ptc.noise[ampName])
                ax.axvline(ptcTurnoff, linestyle="--", color="k", alpha=0.25)
                ax.text(
                    5.5e3, 5.5e4,
                    f"Turnoff: {round(ptcTurnoff, 1)} adu",
                )
                ax.text(
                    5.5e3, 5.1e4,
                    f"Gain: {round(gain, 2)}",
                )
                ax.text(
                    5.5e3, 4.7e4,
                    r"$a_{00}$: " + f"{round(a00*1e6, 2)}" + r"$\times 10^{-6}$",
                )
                ax.text(
                    5.5e3, 4.3e4,
                    r"$n_{00}$: " + f"{round(n00, 2)} el",
                )

                ax.ticklabel_format(style='sci', axis='both', scilimits=(0, 0), useMathText=True)
                ax.set_xlabel("Signal (adu)")
                ax.set_ylabel(r"$C_{00}$ (adu$^2$)")
                # ax.set_xscale('log')#, linthresh=15000.0)
                # ax.set_yscale('log')#, linthresh=15.0)
                ax.set_ylim(0, 6.0e4)
                ax.set_xlim(0, 1.5e5)

            plt.tight_layout()
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}a-"
                        f"ptc_{self.camera[ref.dataId['detector']].getId():03}.png",
                        bbox_inches='tight')
            plt.close()

        for ref in tqdm(self.refs, total=len(self.refs), desc="PTC"):
            ptc = self.butler.get("ptc", dataId=ref.dataId, collections=self.collections)

            fig, axes = plt.subplots(4, 4, figsize=(12, 12))
            fig.suptitle(self.camera[ref.dataId['detector']].getName() + f" ({ref.dataId['detector']})")
            axs = axes.ravel()
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                ampName = amp.getName()
                ax = axs[i]
                ax.set_title(ampName)

                ptcTurnoff = ptc.ptcTurnoff[ampName]

                if ampName in ptc.badAmps:
                    ax.set_facecolor('#fcd4d4')
                    ax.set_title(ampName + " (BAD)")

                data = ptc.finalVars[ampName]
                model = ptc.finalModelVars[ampName]

                ax.scatter(ptc.finalMeans[ampName], data/model - 1, s=1)
                ax.axhline(0, linestyle="--", color="k", alpha=0.25)
                ax.axvline(ptcTurnoff, linestyle="--", color="orange", alpha=0.5)

                ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)
                ax.set_xlabel("Signal (adu)")
                ax.set_ylabel(r"$C_{00} / C_{00}^{model} - 1$")
                ax.set_ylim(-.025, .025)
                ax.set_xlim(0, 1.0e5)
            plt.tight_layout()
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}b-"
                        f"ptc_C00_residuals_{self.camera[ref.dataId['detector']].getId():03}.png",
                        bbox_inches='tight')
            plt.close()

        for ref in tqdm(self.refs, total=len(self.refs), desc="PTC"):
            ptc = self.butler.get("ptc", dataId=ref.dataId, collections=self.collections)

            fig, axes = plt.subplots(4, 4, figsize=(12, 12))
            fig.suptitle(self.camera[ref.dataId['detector']].getName() + f" ({ref.dataId['detector']})")
            axs = axes.ravel()
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                ampName = amp.getName()
                ax = axs[i]
                ax.set_title(ampName)

                ptcTurnoff = ptc.ptcTurnoff[ampName]

                if ampName in ptc.badAmps:
                    ax.set_facecolor('#fcd4d4')
                    ax.set_title(ampName + " (BAD)")

                data = ptc.finalVars[ampName]
                model = ptc.finalModelVars[ampName]

                ax.scatter(ptc.finalMeans[ampName], (data-model)/np.sqrt(data), s=1)
                ax.axvline(ptcTurnoff, linestyle="--", color="orange", alpha=0.5)
                ax.axhline(0, linestyle="--", color="k", alpha=0.25)

                ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)
                ax.set_xlabel("Signal (adu)")
                ax.set_ylabel(r"$C_{00} - C_{00}^{model} / \sqrt{C_{00}}$")
                ax.set_ylim(-3, 3)
                ax.set_xlim(0, 1.0e5)
            plt.tight_layout()
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}c-"
                        f"ptc_C00_residuals_sigma_{self.camera[ref.dataId['detector']].getId():03}.png",
                        bbox_inches='tight')
            plt.close()

        for ref in tqdm(self.refs, total=len(self.refs), desc="PTC"):
            ptc = self.butler.get("ptc", dataId=ref.dataId, collections=self.collections)

            fig, axes = plt.subplots(4, 4, figsize=(12, 12))
            fig.suptitle(self.camera[ref.dataId['detector']].getName() + f" ({ref.dataId['detector']})")
            axs = axes.ravel()
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                ampName = amp.getName()
                ax = axs[i]
                ax.set_title(ampName)

                if ampName in ptc.badAmps:
                    ax.set_facecolor('#fcd4d4')
                    ax.set_title(ampName + " (BAD)")

                ptcTurnoff = ptc.ptcTurnoff[ampName]
                data = ptc.covariances[ampName][:, 1, 0]
                model = ptc.covariancesModel[ampName][:, 1, 0]

                ax.scatter(ptc.finalMeans[ampName], data/model - 1, s=1)
                ax.axhline(0, linestyle="--", color="k", alpha=0.25)
                ax.axvline(ptcTurnoff, linestyle="--", color="orange", alpha=0.5)

                ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)
                ax.set_xlabel("Signal (adu)")
                ax.set_ylabel(r"$C_{10} / C_{10}^{model} - 1$")
                # ax.set_xscale('log')#, linthresh=15000.0)
                # ax.set_yscale('log')#, linthresh=15.0)
                ax.set_ylim(-2, 2)
                ax.set_xlim(0, 1.0e5)
            plt.tight_layout()
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}d-"
                        f"ptc_C10_residuals_{self.camera[ref.dataId['detector']].getId():03}.png",
                        bbox_inches='tight')
            plt.close()

        for ref in tqdm(self.refs, total=len(self.refs), desc="PTC"):
            ptc = self.butler.get("ptc", dataId=ref.dataId, collections=self.collections)

            fig, axes = plt.subplots(4, 4, figsize=(12, 12))
            fig.suptitle(self.camera[ref.dataId['detector']].getName() + f" ({ref.dataId['detector']})")
            axs = axes.ravel()
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                ampName = amp.getName()
                ax = axs[i]
                ax.set_title(ampName)

                if ampName in ptc.badAmps:
                    ax.set_facecolor('#fcd4d4')
                    ax.set_title(ampName + " (BAD)")

                ptcTurnoff = ptc.ptcTurnoff[ampName]
                data = ptc.covariances[ampName][:, 0, 1]
                model = ptc.covariancesModel[ampName][:, 0, 1]

                ax.scatter(ptc.finalMeans[ampName], data/model - 1, s=1)
                ax.axhline(0, linestyle="--", color="k", alpha=0.25)
                ax.axvline(ptcTurnoff, linestyle="--", color="orange", alpha=0.5)

                ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)
                ax.set_xlabel("Signal (adu)")
                ax.set_ylabel(r"$C_{01} / C_{01}^{model} - 1$")
                # ax.set_xscale('log')#, linthresh=15000.0)
                # ax.set_yscale('log')#, linthresh=15.0)
                ax.set_ylim(-1, 1)
                ax.set_xlim(0, 1.0e5)
            plt.tight_layout()
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}e-"
                        f"ptc_C01_residuals_{self.camera[ref.dataId['detector']].getId():03}.png",
                        bbox_inches='tight')
            plt.close()

        for ref in tqdm(self.refs, total=len(self.refs), desc="PTC"):
            ptc = self.butler.get("ptc", dataId=ref.dataId, collections=self.collections)

            fig, axes = plt.subplots(4, 4, figsize=(12, 12))
            fig.suptitle(self.camera[ref.dataId['detector']].getName() + f" ({ref.dataId['detector']})")
            axs = axes.ravel()
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                ampName = amp.getName()
                ax = axs[i]
                ax.set_title(ampName)

                if ampName in ptc.badAmps:
                    ax.set_facecolor('#fcd4d4')
                    ax.set_title(ampName + " (BAD)")

                ptcTurnoff = ptc.ptcTurnoff[ampName]

                data = ptc.rawVars[ampName]/ptc.rawMeans[ampName]
                model = ptc.finalModelVars[ampName] / ptc.finalMeans[ampName]

                ax.scatter(
                    ptc.rawMeans[ampName],
                    ptc.rawVars[ampName]/ptc.rawMeans[ampName],
                    s=1,
                    label="Data",
                )
                ax.plot(ptc.finalMeans[ampName], model, ms=1, color="orange", label="Model")
                ax.axvline(ptcTurnoff, linestyle="--", color="forestgreen", alpha=0.5, label="PTC Turnoff")

                ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)
                ax.set_xlabel("Signal (adu)")
                ax.set_ylabel(r"$C_{00} / \mu$")
                # ax.set_xscale('log')#, linthresh=15000.0)
                # ax.set_yscale('log')#, linthresh=15.0)
                ax.set_ylim(0.4, 0.8)
                ax.set_xlim(0, 1.0e5)
                ax.legend(loc=1, fontsize=8)
            plt.tight_layout()
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}f-"
                        f"ptc_correlations_C00_{self.camera[ref.dataId['detector']].getId():03}.png",
                        bbox_inches='tight')
            plt.close()

        for ref in tqdm(self.refs, total=len(self.refs), desc="PTC"):
            ptc = self.butler.get("ptc", dataId=ref.dataId, collections=self.collections)

            fig, axes = plt.subplots(4, 4, figsize=(12, 12))
            fig.suptitle(self.camera[ref.dataId['detector']].getName() + f" ({ref.dataId['detector']})")
            axs = axes.ravel()
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                ampName = amp.getName()
                ax = axs[i]
                ax.set_title(ampName)

                if ampName in ptc.badAmps:
                    ax.set_facecolor('#fcd4d4')
                    ax.set_title(ampName + " (BAD)")

                ptcTurnoff = ptc.ptcTurnoff[ampName]

                data = ptc.covariances[ampName][:, 1, 0] / ptc.rawMeans[ampName]
                model = ptc.covariancesModel[ampName][:, 1, 0] / ptc.finalMeans[ampName]

                ax.scatter(ptc.rawMeans[ampName], data, s=1, label="Data")
                ax.plot(ptc.finalMeans[ampName], model, ms=1, color="orange", label="Model")
                ax.axhline(0, linestyle="--", color="k", alpha=0.5)
                ax.axvline(ptcTurnoff, linestyle="--", color="forestgreen", alpha=0.5, label="PTC Turnoff")

                ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)
                ax.set_xlabel("Signal (adu)")
                ax.set_ylabel(r"$C_{10} / \mu$")
                # ax.set_xscale('log')#, linthresh=15000.0)
                # ax.set_yscale('log')#, linthresh=15.0)
                ax.set_ylim(-0.01, 0.02)
                ax.set_xlim(0, 1.0e5)
                ax.legend(loc=2, fontsize=8)
            plt.tight_layout()
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}g-"
                        f"ptc_correlations_C10_{self.camera[ref.dataId['detector']].getId():03}.png",
                        bbox_inches='tight')
            plt.close()

        for ref in tqdm(self.refs, total=len(self.refs), desc="PTC"):
            ptc = self.butler.get("ptc", dataId=ref.dataId, collections=self.collections)

            fig, axes = plt.subplots(4, 4, figsize=(12, 12))
            fig.suptitle(self.camera[ref.dataId['detector']].getName() + f" ({ref.dataId['detector']})")
            axs = axes.ravel()
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                ampName = amp.getName()
                ax = axs[i]
                ax.set_title(ampName)

                if ampName in ptc.badAmps:
                    ax.set_facecolor('#fcd4d4')
                    ax.set_title(ampName + " (BAD)")

                ptcTurnoff = ptc.ptcTurnoff[ampName]

                data = ptc.covariances[ampName][:, 0, 1] / ptc.rawMeans[ampName]
                model = ptc.covariancesModel[ampName][:, 0, 1] / ptc.finalMeans[ampName]

                ax.scatter(ptc.rawMeans[ampName], data, s=1, label="Data")
                ax.plot(ptc.finalMeans[ampName], model, ms=1, color="orange", label="Model")
                ax.axhline(0, linestyle="--", color="k", alpha=0.5)
                ax.axvline(ptcTurnoff, linestyle="--", color="forestgreen", alpha=0.5, label="PTC Turnoff")

                ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)
                ax.set_xlabel("Signal (adu)")
                ax.set_ylabel(r"$C_{01} / \mu$")
                # ax.set_xscale('log')#, linthresh=15000.0)
                # ax.set_yscale('log')#, linthresh=15.0)
                ax.set_ylim(-0.01, 0.05)
                ax.set_xlim(0, 1.0e5)
                ax.legend(loc=2, fontsize=8)
            plt.tight_layout()
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}h-"
                        f"ptc_correlations_C01_{self.camera[ref.dataId['detector']].getId():03}.png",
                        bbox_inches='tight')
            plt.close()

        for ref in tqdm(self.refs, total=len(self.refs), desc="PTC"):
            ptc = self.butler.get("ptc", dataId=ref.dataId, collections=self.collections)

            fig, axes = plt.subplots(4, 4, figsize=(12, 12))
            fig.suptitle(self.camera[ref.dataId['detector']].getName() + f" ({ref.dataId['detector']})")
            axs = axes.ravel()
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                ampName = amp.getName()
                ax = axs[i]
                ax.set_title(ampName)

                if ampName in ptc.badAmps:
                    ax.set_facecolor('#fcd4d4')
                    ax.set_title(ampName + " (BAD)")

                ptcTurnoff = ptc.ptcTurnoff[ampName]

                data01 = ptc.covariances[ampName][:, 0, 1]
                model01 = ptc.covariancesModel[ampName][:, 0, 1]
                data10 = ptc.covariances[ampName][:, 0, 1]
                model10 = ptc.covariancesModel[ampName][:, 0, 1]

                ax.scatter(ptc.rawMeans[ampName], data01/data10, s=1, label="Data")
                ax.plot(ptc.finalMeans[ampName], model01/model10, ms=1, color="orange", label="Model")
                ax.axhline(0, linestyle="--", color="k", alpha=0.5)
                ax.axvline(ptcTurnoff, linestyle="--", color="forestgreen", alpha=0.5, label="PTC Turnoff")

                ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)
                ax.set_xlabel("Signal (adu)")
                ax.set_ylabel(r"$C_{01} / C_{10}$")
                # ax.set_xscale('log')#, linthresh=15000.0)
                # ax.set_yscale('log')#, linthresh=15.0)
                ax.set_ylim(-1, 5)
                ax.set_xlim(0, 1.0e5)
                ax.legend(loc=2, fontsize=8)
            plt.tight_layout()
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}i-"
                        f"ptc_C01_C10_ratio_{self.camera[ref.dataId['detector']].getId():03}.png",
                        bbox_inches='tight')
            plt.close()

        for ref in tqdm(self.refs, total=len(self.refs), desc="PTC"):
            ptc = self.butler.get("ptc", dataId=ref.dataId, collections=self.collections)

            fig, axes = plt.subplots(4, 4, figsize=(12, 12))
            fig.suptitle(self.camera[ref.dataId['detector']].getName() + f" ({ref.dataId['detector']})")
            axs = axes.ravel()
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                ampName = amp.getName()
                ax = axs[i]
                ax.set_title(ampName)

                ptcTurnoff = ptc.ptcTurnoff[ampName]

                data = ptc.noiseMatrix[ampName]
                im = ax.imshow(data, origin='lower', norm=SymLogNorm(vmin=0, vmax=15**2, linthresh=1))
                plt.colorbar(im, ax=ax, shrink=0.7, label=r"$n_{ij}$ (electron$^{2}$)")
                ax.set_xlabel("Serial")
                ax.set_ylabel("Parallel")
            plt.tight_layout()
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}j-"
                        f"ptc_noise_matrix_{self.camera[ref.dataId['detector']].getId():03}.png",
                        bbox_inches='tight')
            plt.close()

    def plotBfkStatistics(self):
        self.plotIdx += 1

        # Get datasets, if they exist
        try:
            refs = self.butler.query_datasets(
                "bfk",
                instrument=self.instrument,
                where=f"detector in ({self.detectorIdsString})",
                collections=self.collections,
            )
        except EmptyQueryResultError:
            return

        self.refs = refs

        for ref in tqdm(self.refs, total=len(self.refs), desc="BFK"):
            bfk = self.butler.get("bfk", dataId=ref.dataId, collections=self.collections)
            fig, axes = plt.subplots(4, 4, figsize=(12, 12))
            fig.suptitle(self.camera[ref.dataId['detector']].getName() + f" ({ref.dataId['detector']})")
            axs = axes.ravel()
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                ampName = amp.getName()
                ax = axs[i]

                ax.set_title(ampName)
                if not bfk.valid[ampName]:
                    ax.set_facecolor('#fcd4d4')
                    ax.set_title(ampName + " (INVALID)")

                length = bfk.ampKernels[ampName].shape[0]
                ax.plot(bfk.ampKernels[ampName][length//2, :], drawstyle="steps-mid", label="Serial\nslice")
                ax.plot(bfk.ampKernels[ampName][:, length//2], drawstyle="steps-mid", label="Parallel\nslice")
                ax.set_xticks([0, length//2, length-1], [-length//2 + 1, 0, length//2])

                ax.axhline(0, linestyle="--", color="k", alpha=0.5)

                ax.set_xlabel("Offset (px)")
                ax.set_ylabel(r"(unitless)")
                ax.set_ylim(-2.5e-6, 0.25e-6)
                ax.legend(loc=3, fontsize=8)
            plt.tight_layout()
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}-"
                        f"bfk_profiles_{self.camera[ref.dataId['detector']].getId():03}.png",
                        bbox_inches='tight')
            plt.close()

    def plotCtiStatistics(self):
        self.plotIdx += 1

        # Get datasets, if they exist
        try:
            refs = self.butler.query_datasets(
                "cti",
                instrument=self.instrument,
                where=f"detector in ({self.detectorIdsString})",
                collections=self.collections,
            )
        except EmptyQueryResultError:
            return

        self.refs = refs

        for ref in tqdm(self.refs, total=len(self.refs), desc="CTI"):
            cti = self.butler.get("cti", dataId=ref.dataId, collections=self.collections)
            lin = self.butler.get("linearizer", dataId=ref.dataId, collections=self.collections)
            ptc = self.butler.get("ptc", dataId=ref.dataId, collections=self.collections)

            fig, axes = plt.subplots(4, 4, figsize=(12, 12))
            fig.suptitle(self.camera[ref.dataId['detector']].getName())
            axs = axes.ravel()
            detName = self.camera[ref.dataId['detector'].getName()]
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                ampName = amp.getName()

                # Save dataset information
                self.typesDict[detName][ampName] = self.camera[ref.dataId['detector']].getPhysicalType()
                self.linearityTurnoffsDict[detName][ampName] = lin.linearityTurnoff[ampName]
                self.ptcTurnoffsDict[detName][ampName] = ptc.ptcTurnoff[ampName] * ptc.gain[ampName]
                self.globalCtisDict[detName][ampName] = cti.globalCti[ampName]
                self.serialCtiTurnoffsDict[detName][ampName] = cti.serialCtiTurnoff[ampName]

                ax = axs[i]
                ax.set_title(ampName)
                ax.scatter(cti.signals[ampName], cti.serialEper[ampName], s=1)

                globalCti = cti.globalCti[ampName]
                serialCtiTurnoff = cti.serialCtiTurnoff[ampName]

                ax.axhline(0, linestyle="--", color="k", alpha=0.5)

                try:
                    power = int(np.floor(np.log10(np.abs(globalCti))))
                    num = globalCti * (10**np.abs(power))
                except OverflowError:
                    power = 0
                    num = 0

                ax.axhline(globalCti, linestyle="--", color="orange", alpha=0.5,
                           label=f"Global sCTI: {round(num, 3)}" + rf"$\times10^{{{power}}}$")
                ax.axhline(5e-6, linestyle="-", color="red", alpha=0.25, label="Max spec.")
                ax.axvline(serialCtiTurnoff, linestyle="--", color="forestgreen", alpha=0.5,
                           label=f"sCTI Turnoff: {round(serialCtiTurnoff, 1)}")

                ax.ticklabel_format(style='sci', axis='both', scilimits=(0, 0), useMathText=True)
                ax.set_xlabel("Signal (el)")
                ax.set_ylabel(r"Serial EPER")
                ax.set_yscale('log')
                ax.set_ylim(1e-7, 1e-5)
                ax.set_xlim(0, 2.0e5)
                ax.legend(fontsize=8)
            plt.tight_layout()
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}-"
                        f"serial_eper_{self.camera[ref.dataId['detector']].getId():03}.png",
                        bbox_inches='tight')
            plt.close()

    def plotBiasStatistics(self):
        self.plotIdx += 1

        # Get dataset references, if they exist
        try:
            refs = self.butler.query_datasets(
                "bias",
                instrument=self.instrument,
                where=f"detector in ({self.detectorIdsString})",
                collections=self.collections,
            )
        except EmptyQueryResultError:
            return

        self.refs = refs

        for ref in tqdm(self.refs, total=len(self.refs), desc="BIAS"):
            plt.figure()
            img = self.butler.get("bias", dataId=ref.dataId, collections=self.collections)
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                bbox = amp.getBBox()
                vals = img.image[bbox].array.ravel()
                vals = vals[(vals > -10.0) * (vals < 10.0)]
                _ = plt.hist(vals, bins="fd", histtype='step', color=self.color[i], label=amp.getName())

            plt.yscale('log')
            plt.xlabel("electrons")
            plt.legend(loc=2, ncol=2)
            plt.xlim(-10, 10)
            plt.axvline(0.0, linestyle="--", color="k", alpha=0.5)
            plt.title(f"BIAS {self.camera[ref.dataId['detector']].getName()}, {ref.dataId['detector']}")
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}-"
                        f"bias_histogram_{ref.dataId['detector']:03}.png",
                        bbox_inches='tight')
            plt.close()

    def plotDarkStatistics(self):
        self.plotIdx += 1

        # Get dataset references, if they exist
        try:
            refs = self.butler.query_datasets(
                "dark",
                instrument=self.instrument,
                where=f"detector in ({self.detectorIdsString})",
                collections=self.collections,
            )
        except EmptyQueryResultError:
            return

        self.refs = refs

        for ref in tqdm(self.refs, total=len(self.refs), desc="DARK"):
            plt.figure()
            img = self.butler.get("dark", dataId=ref.dataId, collections=self.collections)
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                bbox = amp.getBBox()
                vals = img.image[bbox].array.ravel()
                vals = vals[(vals > -2.5) * (vals < 2.5)]
                _ = plt.hist(vals, bins="fd", histtype='step', color=self.color[i], label=amp.getName())

            plt.yscale('log')
            plt.xlabel("electrons")
            plt.legend(loc=2, ncol=2)
            plt.xlim(-2.5, 2.5)
            plt.axvline(0.0, linestyle="--", color="k", alpha=0.5)
            plt.title(f"DARK {self.camera[ref.dataId['detector']].getName()}, {ref.dataId['detector']}")
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}-"
                        f"dark_histogram_{ref.dataId['detector']:03}.png",
                        bbox_inches='tight')
            plt.close()

    def plotBiasVsDarkStatistics(self):
        self.plotIdx += 1

        if self.refs is None:
            self.log.warning("There are no bias and darks, skipping plotBiasVsDarkStatistics.")

        for ref in tqdm(self.refs, total=len(self.refs), desc="BIAS VS DARK"):
            try:
                img = self.butler.get(
                    "bias",
                    dataId=ref.dataId,
                    collections=self.collections
                )
                img2 = self.butler.get(
                    "dark",
                    dataId=ref.dataId,
                    collections=self.collections
                )
            except DatasetNotFoundError:
                self.log.warning(f"No matching bias and dark datasets for {ref}, skipping.")
                continue

            fig, axes = plt.subplots(4, 4, figsize=(10, 10))
            fig.suptitle(self.camera[ref.dataId['detector']].getName() + f" ({ref.dataId['detector']})")
            axs = axes.ravel()
            for i, amp in enumerate(self.camera[ref.dataId['detector']].getAmplifiers()):
                ax = axs[i]
                bbox = amp.getBBox()
                vals = img.image[bbox].array.ravel()
                vals2 = img2.image[bbox].array.ravel() * 30.0  # scaled by 30s
                ax.hexbin(vals, vals2, bins='log', gridsize=(100, 100), extent=(-15, 15, -30, 30))
                ax.set_xlabel("BIAS")
                ax.set_ylabel("DARK * 30s")
                ax.set_xlim(-15., 15.)
                ax.set_ylim(-30., 30.)
                ax.set_title(amp.getName())
                ax.set_aspect(0.5)
            plt.tight_layout()
            outputFigPath = f"{self.output}/{self._getRaftBayPathFromRef(ref)}"
            plt.savefig(f"{outputFigPath}/step_{self.plotIdx:02}-"
                        f"bias_and_dark_{ref.dataId['detector']:03}.png",
                        bbox_inches='tight')
            plt.close()

    def plotSummaryStatistics(self):
        if len(self.refs) > 0:
            self.log.info("Plotting summary statistics")
        else:
            self.log.warning("No data for summary statistics, skipping.")
            return

        # Reset plot index to zero
        self.plotIdx = 0

        # Get the physical type masks
        types = self._getAllDictValues(self.typesDict)
        e2vs = (types == "E2V")
        itls = (types == "ITL")

        self.plotIdx += 1
        if np.any(self._getAllDictValues(self.n00sDict)):
            plt.figure(figsize=(12, 12))
            plot_fp_statistic(plt.gca(), self.n00sDict, camera=self.camera, cm=plt.cm.hot,
                              x_range=None, y_range=None,
                              z_range=(4, 15), use_log10=False, scale_factor='1',
                              title=r'$n_{00}$ (electrons)', nsigma=4)
            plt.savefig(f"{self.output}/{self.plotIdx:02}a-n00_mosaic.png",
                        bbox_inches='tight')
            plt.close()

            plt.figure()
            n00s = self._getAllDictValues(self.n00sDict)
            plt.hist(n00s[e2vs], bins="fd", histtype="step", label="E2Vs")
            plt.hist(n00s[itls], bins="fd", histtype="step", label="ITLs")
            plt.legend(loc=1)
            plt.yscale('log')
            plt.xlabel("noise (electron)")
            plt.savefig(f"{self.output}/{self.plotIdx:02}b-n00_hist_mosaic.png",
                        bbox_inches='tight')
            plt.close()

        self.plotIdx += 1
        if np.any(self._getAllDictValues(self.empiricalN00sDict)):
            plt.figure(figsize=(12, 12))

            plot_fp_statistic(plt.gca(), self.empiricalN00sDict, camera=self.camera, cm=plt.cm.hot,
                              x_range=None, y_range=None,
                              z_range=(4, 15), use_log10=False, scale_factor='1',
                              title='Mean Overscan Residual Sigma (Input PTC Flats)', nsigma=4)
            plt.savefig(f"{self.output}/{self.plotIdx:02}a-overscan_residuals_sigma_mosaic.png",
                        bbox_inches='tight')
            plt.close()

            plt.figure()
            empiricalN00s = self._getAllDictValues(self.empiricalN00sDict)
            plt.hist(empiricalN00s[e2vs], bins="fd", histtype="step", label="E2Vs")
            plt.hist(empiricalN00s[itls], bins="fd", histtype="step", label="ITLs")
            plt.legend(loc=1)
            plt.yscale('log')
            plt.xlabel("overscan residual sigma\n(median of used input PTC flats)", fontsize=12)
            plt.savefig(f"{self.output}/{self.plotIdx:02}b-overscan_residual_sigma_hist_mosaic.png",
                        bbox_inches='tight')
            plt.close()

        self.plotIdx += 1
        if np.any(self._getAllDictValues(self.empiricalN00sResidualsDict)):
            plt.figure(figsize=(12, 12))
            plot_fp_statistic(plt.gca(), self.empiricalN00sResidualsDict, camera=self.camera, cm=plt.cm.bwr,
                              x_range=None, y_range=None,
                              z_range=(-0.15, 0.15), use_log10=False, scale_factor='1',
                              title='Noise / Mean Overscan Residual Sigma - 1\n(Input PTC Flats)', nsigma=4)
            plt.savefig(f"{self.output}/{self.plotIdx:02}a-noise_and_overscan_residuals_sigma_mosaic.png",
                        bbox_inches='tight')
            plt.close()

            plt.figure()
            empiricalN00sResiduals = self._getAllDictValues(self.empiricalN00sResidualsDict)
            plt.hist(empiricalN00sResiduals[e2vs], bins=100, histtype="step", label="E2Vs")
            plt.hist(empiricalN00sResiduals[itls], bins=100, histtype="step", label="ITLs")
            plt.yscale('log')
            plt.legend(loc=1)
            plt.xlabel("noise / oscanResidStd - 1", fontsize=12)
            name = "noise_and_overscan_residual_sigma_residuals_hist_mosaic"
            plt.savefig(f"{self.output}/{self.plotIdx:02}b-{name}.png",
                        bbox_inches='tight')
            plt.close()

        self.plotIdx += 1
        if np.any(self._getAllDictValues(self.gainsDict)):
            plt.figure(figsize=(12, 12))
            plot_fp_statistic(plt.gca(), self.gainsDict, camera=self.camera, cm=plt.cm.hot,
                              x_range=None, y_range=None,
                              z_range="clipped_autoscale", use_log10=False, scale_factor='1',
                              title='Gains', nsigma=4)
            plt.savefig(f"{self.output}/{self.plotIdx:02}a-gains_mosaic.png",
                        bbox_inches='tight')
            plt.close()

            plt.figure()
            gains = self._getAllDictValues(self.gainsDict)
            plt.hist(gains[e2vs], bins="fd", histtype="step", label="E2Vs")
            plt.hist(gains[itls], bins="fd", histtype="step", label="ITLs")
            plt.yscale('log')
            plt.legend(loc=1)
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)
            plt.xlabel(r"Gains (el/adu)")
            plt.savefig(f"{self.output}/{self.plotIdx:02}b-gains_hist_mosaic.png",
                        bbox_inches='tight')
            plt.close()

        self.plotIdx += 1
        if np.any(self._getAllDictValues(self.gainsUnadjustedDict)):
            plt.figure(figsize=(12, 12))
            plot_fp_statistic(plt.gca(), self.gainsUnadjustedDict, camera=self.camera, cm=plt.cm.hot,
                              x_range=None, y_range=None,
                              z_range="clipped_autoscale", use_log10=False, scale_factor='1',
                              title='Gains (Unadjusted)', nsigma=4)
            plt.savefig(f"{self.output}/{self.plotIdx:02}a-gains_unadjusted_mosaic.png",
                        bbox_inches='tight')
            plt.close()

            plt.figure()
            gainsUnadjusted = self._getAllDictValues(self.gainsUnadjustedDict)
            plt.hist(gainsUnadjusted[e2vs], bins="fd", histtype="step", label="E2Vs")
            plt.hist(gainsUnadjusted[itls], bins="fd", histtype="step", label="ITLs")
            plt.legend(loc=1)
            plt.xlabel(r"Gains (Unadjusted)")
            plt.savefig(f"{self.output}/{self.plotIdx:02}b-gains_unadjusted_hist_mosaic.png",
                        bbox_inches='tight')
            plt.close()

        self.plotIdx += 1
        if np.any(self._getAllDictValues(self.gainsUnadjustedResidualsDict)):
            plt.figure(figsize=(12, 12))
            plot_fp_statistic(plt.gca(), self.gainsUnadjustedResidualsDict, camera=self.camera, cm=plt.cm.bwr,
                              x_range=None, y_range=None,
                              z_range=(-0.010, 0.010), use_log10=False, scale_factor='1',
                              title='Gains (Unadjusted) / Gains - 1', nsigma=4)
            plt.savefig(f"{self.output}/{self.plotIdx:02}a-gains_unadjusted_residuals_mosaic.png",
                        bbox_inches='tight')
            plt.close()

            plt.figure()
            gainsUnadjustedResiduals = self._getAllDictValues(self.gainsUnadjustedResidualsDict)
            plt.hist(gainsUnadjustedResiduals[e2vs], bins="fd", histtype="step", label="E2Vs")
            plt.hist(gainsUnadjustedResiduals[itls], bins="fd", histtype="step", label="ITLs")
            plt.yscale('log')
            plt.legend(loc=1)
            plt.xlabel(r"Gains Unadjusted / Gains - 1")
            plt.axvline(0, linestyle="--", color="k", alpha=0.25)
            plt.savefig(f"{self.output}/{self.plotIdx:02}b-gains_unadjusted_residuals_hist_mosaic.png",
                        bbox_inches='tight')
            plt.close()

        self.plotIdx += 1
        if np.any(self._getAllDictValues(self.a00sDict)):
            plt.figure(figsize=(12, 12))
            plot_fp_statistic(plt.gca(), self.a00sDict, camera=self.camera, cm=plt.cm.hot,
                              x_range=None, y_range=None,
                              z_range="clipped_autoscale", use_log10=False, scale_factor='1',
                              title=r'$a_{00}$', nsigma=4)
            plt.savefig(f"{self.output}/{self.plotIdx:02}a-a00_mosaic.png",
                        bbox_inches='tight')
            plt.close()

            plt.figure()
            a00s = self._getAllDictValues(self.a00sDict)
            plt.hist(a00s[e2vs], bins="fd", histtype="step", label="E2Vs")
            plt.hist(a00s[itls], bins="fd", histtype="step", label="ITLs")
            plt.yscale('log')
            plt.legend(loc=1)
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)
            plt.xlabel(r"$a_{00}$")
            plt.savefig(f"{self.output}/{self.plotIdx:02}b-a00_hist_mosaic.png",
                        bbox_inches='tight')
            plt.close()

        self.plotIdx += 1
        if np.any(self._getAllDictValues(self.n01sDict)):
            plt.figure(figsize=(12, 12))
            plot_fp_statistic(plt.gca(), self.n01sDict, camera=self.camera, cm=plt.cm.hot,
                              x_range=None, y_range=None,
                              z_range=(0, 10), use_log10=False, scale_factor='1',
                              title=r'$n_{01}$ (electrons$^{2}$)', nsigma=5)
            plt.savefig(f"{self.output}/{self.plotIdx:02}a-n01_mosaic.png",
                        bbox_inches='tight')
            plt.close()

            plt.figure()
            n01s = self._getAllDictValues(self.n01sDict)
            plt.hist(n01s[e2vs], bins="fd", histtype="step", label=r"E2Vs")
            plt.hist(n01s[itls], bins="fd", histtype="step", label=r"ITLs")
            plt.yscale('log')
            plt.legend(loc=1)
            plt.xlabel(r"$n_{01}$ (electrons$^2$)")
            plt.savefig(f"{self.output}/{self.plotIdx:02}b-n01_hist_mosaic.png",
                        bbox_inches='tight')
            plt.close()

        self.plotIdx += 1
        if np.any(self._getAllDictValues(self.n10sDict)):
            plt.figure(figsize=(12, 12))
            plot_fp_statistic(plt.gca(), self.n10sDict, camera=self.camera, cm=plt.cm.hot,
                              x_range=None, y_range=None,
                              z_range=(-5, 60), use_log10=False, scale_factor='1',
                              title=r'$n_{10}$ (electrons$^{2}$)', nsigma=4)
            plt.savefig(f"{self.output}/{self.plotIdx:02}a-n10_mosaic.png",
                        bbox_inches='tight')
            plt.close()

            plt.figure()
            n10s = self._getAllDictValues(self.n10sDict)
            plt.hist(n10s[e2vs], bins="fd", histtype="step", label=r"E2Vs")
            plt.hist(n10s[itls], bins="fd", histtype="step", label=r"ITLs")
            plt.yscale('log')
            plt.legend(loc=2)
            plt.xlabel(r"$n_{10}$ (electrons$^2$)")
            plt.savefig(f"{self.output}/{self.plotIdx:02}b-n10_hist_mosaic.png",
                        bbox_inches='tight')
            plt.close()

        self.plotIdx += 1
        if np.any(self._getAllDictValues(self.ptcTurnoffsDict)):
            plt.figure(figsize=(12, 12))
            plot_fp_statistic(plt.gca(), self.ptcTurnoffsDict, camera=self.camera, cm=plt.cm.hot,
                              x_range=None, y_range=None,
                              z_range="clipped_autoscale", use_log10=False, scale_factor='1',
                              title='PTC Turnoffs (electrons)', nsigma=4)
            plt.savefig(f"{self.output}/{self.plotIdx:02}a-ptc_turnoffs_mosaic.png",
                        bbox_inches='tight')
            plt.close()

            plt.figure()
            ptcTurnoffs = self._getAllDictValues(self.ptcTurnoffsDict)
            plt.hist(ptcTurnoffs[e2vs], bins="fd", histtype="step", label="E2Vs")
            plt.hist(ptcTurnoffs[itls], bins="fd", histtype="step", label="ITLs")
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)
            plt.xlabel("PTC Turnoffs (electrons)")
            plt.yscale('log')
            plt.legend(loc=2)
            plt.savefig(f"{self.output}/{self.plotIdx:02}b-ptc_turnoffs_hist_mosaic.png",
                        bbox_inches='tight')
            plt.close()

        self.plotIdx += 1
        if np.any(self._getAllDictValues(self.linearityTurnoffsDict)):
            plt.figure(figsize=(12, 12))
            plot_fp_statistic(plt.gca(), self.linearityTurnoffsDict, camera=self.camera, cm=plt.cm.hot,
                              x_range=None, y_range=None,
                              z_range=None, use_log10=False, scale_factor='1',
                              title=r'Linearity Turnoffs (adu)', nsigma=4)
            plt.savefig(f"{self.output}/{self.plotIdx:02}a-linearity_turnoffs_mosaic.png",
                        bbox_inches='tight')
            plt.close()

            plt.figure()
            linearityTurnoffs = self._getAllDictValues(self.linearityTurnoffsDict)
            plt.hist(linearityTurnoffs[e2vs], bins="fd", histtype="step", label="E2Vs")
            plt.hist(linearityTurnoffs[itls], bins="fd", histtype="step", label="ITLs")
            plt.yscale('log')
            plt.legend(loc=2)
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)
            plt.xlabel("Linearity Turnoffs (adu)", fontsize=12)
            plt.savefig(f"{self.output}/{self.plotIdx:02}b-linearity_turnoffs_hist_mosaic.png",
                        bbox_inches='tight')
            plt.close()

        self.plotIdx += 1
        if np.any(self._getAllDictValues(self.empiricalN00sResidualsDict)):
            plt.figure(figsize=(12, 12))
            plot_fp_statistic(plt.gca(), self.globalCtisDict, camera=self.camera, cm=plt.cm.hot,
                              x_range=None, y_range=None,
                              z_range=(1e-7, 1e-5), use_log10=False, scale_factor='1',
                              title=r'Global sCTI ($transfers^{-1}$)', nsigma=4)
            plt.savefig(f"{self.output}/{self.plotIdx:02}a-global_cti_mosaic.png", bbox_inches='tight')
            plt.close()

            plt.figure()
            globalCtis = self._getAllDictValues(self.globalCtisDict)
            plt.hist(globalCtis[e2vs], bins=100, histtype='step', range=(1e-7, 1e-5), label="E2V")
            plt.hist(globalCtis[itls], bins=100, histtype='step', range=(1e-7, 1e-5), label="ITL")
            plt.axvline(5e-6, linestyle="-", color="r", alpha=0.25, label="Max. Spec.")
            bad = np.count_nonzero(globalCtis >= 5e-6)
            plt.hist([], bins=1, color="white", histtype='step',
                     label=f"({bad}/{len(globalCtis)} out of spec.)")
            plt.xlim(1e-7, 1e-5)
            plt.legend(loc=1)
            plt.yscale('log')
            plt.legend(loc=1)
            plt.xlabel(r"Global sCTI [transfers$^{-1}$]")
            plt.savefig(f"{self.output}/{self.plotIdx:02}b-global_cti_hist_mosaic.png",
                        bbox_inches='tight')
            plt.close()

        self.plotIdx += 1
        if np.any(self._getAllDictValues(self.serialCtiTurnoffsDict)):
            plt.figure(figsize=(12, 12))
            plot_fp_statistic(plt.gca(), self.serialCtiTurnoffsDict, camera=self.camera, cm=plt.cm.hot,
                              x_range=None, y_range=None,
                              z_range="clipped_autoscale", use_log10=False, scale_factor='1',
                              title='Serial CTI Turnoffs (electrons)', nsigma=4)
            plt.savefig(f"{self.output}/{self.plotIdx:02}a-serial_cti_turnoffs_mosaic.png",
                        bbox_inches='tight')
            plt.close()

            plt.figure()
            serialCtiTurnoffs = self._getAllDictValues(self.serialCtiTurnoffsDict)
            plt.hist(serialCtiTurnoffs[e2vs], bins="fd", histtype='step', label="E2V")
            plt.hist(serialCtiTurnoffs[itls], bins="fd", histtype='step', label="ITL")
            plt.legend(loc=2)
            plt.yscale('log')
            plt.xlabel("Serial CTI Turnoffs (electrons)")
            plt.savefig(f"{self.output}/{self.plotIdx:02}b-serial_cti_turnoffs_hist_mosaic.png",
                        bbox_inches='tight')
            plt.close()

    def run(self):
        # Setup
        self.initializeOutputDirectories()

        # Plot
        with logging_redirect_tqdm():
            self.plotLinearizerStatistics()
            self.plotPtcStatistics()
            self.plotBfkStatistics()
            self.plotCtiStatistics()
            self.plotBiasStatistics()
            self.plotDarkStatistics()
            self.plotBiasVsDarkStatistics()
            if self.plotSummaryStatistics:
                self.plotSummaryStatistics()


def main():
    # Could these use the click/daf_butler CLI tools?
    parser = argparse.ArgumentParser(description="Construct a detailed report.")
    parser.add_argument("-r", "--repository", dest="repository", default="",
                        help="Butler repository to pull results from.")
    parser.add_argument("-o", "--output_path", dest="output_path", default="",
                        help="Output path to write report to.")
    parser.add_argument("-c", "--collections", dest="collections", action="append", default=[],
                        help="Collections to search for results.")
    parser.add_argument("--instrument", dest="instrument", default="",
                        help="The instrument to use (e.g. LSSTCam)")
    parser.add_argument("--detector-ids", nargs="+", default=range(205), type=int,
                        help="Specific detectors to plot, default is all available.")
    parser.add_argument("--plot-summary-stats", dest="plot_summary_stats", action="store_true",
                        default=True, help="Plot the summary statistics?")
    parser.add_argument("--do-overwrite", dest="do_overwrite", action="store_true", default=False,
                        help="Allow existing files to be overwritten?")
    args = parser.parse_args()

    plotter = DetailedPlotter(
        repo=args.repository,
        output=args.output_path,
        collections=args.collections,
        instrument=args.instrument,
        detectorIds=args.detector_ids,
        plotSummaryStats=args.plot_summary_stats,
        doOverwrite=args.do_overwrite,
    )
    plotter.run()
