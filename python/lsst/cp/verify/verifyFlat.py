# This file is part of cp_verify.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
import numpy as np

from astropy.table import Table

import lsst.afw.math as afwMath
from .verifyStats import CpVerifyStatsConfig, CpVerifyStatsTask, CpVerifyStatsConnections
from .mergeResults import CpVerifyExpMergeConfig, CpVerifyExpMergeTask

__all__ = ['CpVerifyFlatConfig', 'CpVerifyFlatTask',
           'CpVerifyFlatExpMergeConfig', 'CpVerifyFlatExpMergeTask']
# 'CpVerifyFlatRunMergeConfig', 'CpVerifyFlatRunMergeTask',]


class CpVerifyFlatConfig(CpVerifyStatsConfig,
                         pipelineConnections=CpVerifyStatsConnections):
    """Inherits from base CpVerifyStatsConfig.
    """

    def setDefaults(self):
        super().setDefaults()
        self.imageStatKeywords = {'MEAN': 'MEAN',  # noqa F841
                                  'NOISE': 'STDEVCLIP', }
        self.detectorStatKeywords = {'MEAN': 'MEAN',  # noqa F841
                                     'SCATTER': 'STDEV', }


class CpVerifyFlatTask(CpVerifyStatsTask):
    """Flat verification sub-class, implementing the verify method.
    """
    ConfigClass = CpVerifyFlatConfig
    _DefaultName = 'cpVerifyFlat'

    stageName = "flat"

    def detectorStatistics(self, statisticsDict, statControl, exposure=None, uncorrectedExposure=None):
        """Calculate detector level statistics based on the existing
        per-amplifier measurements.

        Parameters
        ----------
        statisticsDict : `dict` [`str`, `dict` [`str`, scalar]],
            Dictionary of measured statistics.  The inner dictionary
            should have keys that are statistic names (`str`) with
            values that are some sort of scalar (`int` or `float` are
            the mostly likely types).
        statControl : `lsst.afw.math.StatControl`
            Statistics control object with parameters defined by
            the config.
        exposure : `lsst.afw.image.Exposure`, optional
            Exposure containing the ISR-processed data to measure.
        uncorrectedExposure : `lsst.afw.image.Exposure`, optional
            uncorrected esposure (no defects) containing the
            ISR-processed data to measure.

        Returns
        -------
        outputStatistics : `dict` [`str`, scalar]
            A dictionary of the statistics measured and their values.

        """
        outputStatistics = {}

        ampStats = statisticsDict['AMP']
        amplifierMeans = [stats['MEAN'] for stats in ampStats.values()]

        statisticToRun, statAccessor = self._configHelper(self.config.detectorStatKeywords)
        stats = afwMath.makeStatistics(amplifierMeans, statisticToRun, statControl)

        for k, v in statAccessor.items():
            outputStatistics[k] = stats.getValue(v)

        return outputStatistics

    def verify(self, exposure, statisticsDict):
        """Verify that the measured statistics meet the verification criteria.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure the statistics are from.
        statisticsDict : `dict` [`str`, `dict` [`str`, scalar]],
            Dictionary of measured statistics.  The inner dictionary
            should have keys that are statistic names (`str`) with
            values that are some sort of scalar (`int` or `float` are
            the mostly likely types).

        Returns
        -------
        outputStatistics : `dict` [`str`, `dict` [`str`, `bool`]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the verification criteria.
        success : `bool`
            A boolean indicating if all tests have passed.
        """
        ampStats = statisticsDict['AMP']
        verifyStats = {}
        success = True
        for ampName, stats in ampStats.items():
            verify = {}

            # DMTN-101 Test 10.X: confirm that per-amplifier scatter is
            #                     consistent with Poissonian
            verify['NOISE'] = bool(stats['NOISE'] <= np.sqrt(stats['MEAN']))

            verify['SUCCESS'] = bool(np.all(list(verify.values())))
            if verify['SUCCESS'] is False:
                success = False

            verifyStats[ampName] = verify

        verifyDet = {}
        detStats = statisticsDict['DET']

        # DMTN-101 Test 10.Y: confirm intra-chip scatter is small.
        verifyDet['SCATTER'] = bool(detStats['SCATTER']/detStats['MEAN'] <= 0.05)

        verifyDet['SUCCESS'] = bool(np.all(list(verifyDet.values())))
        if verifyDet['SUCCESS'] is False:
            success = False

        return {'AMP': verifyStats, 'DET': verifyDet}, bool(success)

    def repackStats(self, detStats, detDims):
        """Repack hierarchical results into flat table.

        Parameters
        ----------
        detStats : `dict` [`str`, `dict`]
            A nested set of dictionaries containing relevant
            statistics.
        detDims : `dict` [`str`, `str`]
            The dimensions of this set of statistics.

        Returns
        -------
        outputResults : `astropy.Table`, optional
            The repacked flat table.
        outputMatrix : `astropy.Table`, optional
            The repackaed matrix data, in a flat table.
        """
        rowList = []

        row = {}
        instrument = detDims["instrument"]
        exposure = detDims["exposure"]
        detector = detDims["detector"]
        mjd = detStats["ISR"]["MJD"]

        # Get amp stats
        # AMP {ampName} [MEAN NOISE] value
        for ampName, stats in detStats["AMP"].items():
            row[ampName] = {
                "instrument": instrument,
                "exposure": exposure,
                "detector": detector,
                "amplifier": ampName,
                "mjd": mjd,
                "flatMean": stats["MEAN"],
                "flatNoise": stats["NOISE"],
            }
        # Get catalog stats CATALOG

        # Get metadata stats
        # METADATA (RESIDUAL STDEV) {ampName} value

        # Get verify stats
        for ampName, stats in detStats["VERIFY"]["AMP"].items():
            row[ampName]["flatVerifyNoise"] = stats["NOISE"]

        # Get isr stats
        for ampName, stats in detStats["ISR"]["CALIBDIST"].items():
            for level in self.config.expectedDistributionLevels:
                key = f"LSST CALIB {self.stageName.upper()} {ampName} DISTRIBUTION {level}-PCT"
                row[ampName][f"flatDistribution_{level}"] = stats[key]
        # Get detector stats
        # DET
        row["detector"] = {
            "instrument": instrument,
            "exposure": exposure,
            "detector": detector,
            "flatDetMean": detStats["DET"]["MEAN"],
            "flatDetScatter": detStats["DET"]["SCATTER"],
        }

        # Append to output
        for ampName, stats in row.items():
            rowList.append(stats)

        return Table(rowList), None


class CpVerifyFlatExpMergeConfig(CpVerifyExpMergeConfig):
    """Inherits from base CpVerifyExpMergeConfig
    """

    def setDefaults(self):
        super().setDefaults()
        self.exposureStatKeywords = {'EXPOSURE_SCATTER': 'STDEV',  # noqa F841
                                     }


class CpVerifyFlatExpMergeTask(CpVerifyExpMergeTask):
    """Inherits from base CpVerifyExpMergeTask
    """

    def exposureStatistics(self, statisticsDictionary):
        """Calculate exposure level statistics based on the existing
        per-amplifier and per-detector measurements.

        Parameters
        ----------
        statisticsDictionary : `dict [`str`, `dict` [`str`, scalar]],
            Dictionary of measured statistics.  The top level
            dictionary is keyed on the detector names, and contains
            the measured statistics from the per-detector
            measurements.

        Returns
        -------
        outputStatistics : `dict` [`str, scalar]
            A dictionary of the statistics measured and their values.
        """
        detectorMeans = []
        for detName, stats in statisticsDictionary.items():
            # Get detector stats:
            detectorMeans.append(stats['DET']['MEAN'])

        return {'SCATTER': np.stdev(detectorMeans)}

    def verify(self, detectorStatistics, statisticsDictionary):
        """Verify if the measured statistics meet the verification criteria.

        Parameters
        ----------
        detectorStatistics : `dict` [`str`, `dict` [`str`, scalar]]
            Merged set of input detector level statistics.
        statisticsDictionary : `dict` [`str`, `dict` [`str`, scalar]],
            Dictionary of measured statistics.  The inner dictionary
            should have keys that are statistic names (`str`) with
            values that are some sort of scalar (`int` or `float` are
            the mostly likely types).

        Returns
        -------
        outputStatistics : `dict` [`str`, `dict` [`str`, `bool`]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the verification criteria.
        success : `bool`
            A boolean indicating if all tests have passed.

        """
        verifyStats = {}
        success = True

        # DMTN-101 Test 10.Z: confirm inter-chip scatter is small.
        verifyStats['SCATTER'] = bool(statisticsDictionary['EXP']['SCATTER'] <= 0.05)

        success = bool(np.all(list(verifyStats.values())))

        return {'EXP': verifyStats}, bool(success)
