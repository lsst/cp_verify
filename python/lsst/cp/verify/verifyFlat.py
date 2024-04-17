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
import lsst.afw.math as afwMath
from .verifyStats import CpVerifyStatsConfig, CpVerifyStatsTask, CpVerifyStatsConnections
from .mergeResults import CpVerifyExpMergeByFilterConfig, CpVerifyExpMergeByFilterTask

__all__ = ['CpVerifyFlatConfig', 'CpVerifyFlatTask',
           'CpVerifyFlatExpMergeConfig', 'CpVerifyFlatExpMergeTask']
# 'CpVerifyFlatRunMergeConfig', 'CpVerifyFlatRunMergeTask',]


class CpVerifyFlatConfig(CpVerifyStatsConfig,
                         pipelineConnections=CpVerifyStatsConnections):
    """Inherits from base CpVerifyStatsConfig.
    """

    def setDefaults(self):
        super().setDefaults()
        self.stageName = 'FLAT'
        self.imageStatKeywords = {'MEAN': 'MEAN',  # noqa F841
                                  'NOISE': 'STDEVCLIP', }
        self.detectorStatKeywords = {'MEAN': 'MEAN',  # noqa F841
                                     'SCATTER': 'STDEV', }


class CpVerifyFlatTask(CpVerifyStatsTask):
    """Flat verification sub-class, implementing the verify method.
    """
    ConfigClass = CpVerifyFlatConfig
    _DefaultName = 'cpVerifyFlat'

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

    def repackStats(self, statisticsDict, dimensions):
        # docstring inherited
        rows = {}
        rowList = []
        matrixRowList = None

        if self.config.useIsrStatistics:
            mjd = statisticsDict["ISR"]["MJD"]
        else:
            mjd = np.nan

        print(dimensions)
        rowBase = {
            "instrument": dimensions["instrument"],
            "exposure": dimensions["exposure"],
            "detector": dimensions["detector"],
            "physical_filter": dimensions["physical_filter"],
            "mjd": mjd,
        }

        # AMP results:
        for ampName, stats in statisticsDict["AMP"].items():
            rows[ampName] = {}
            rows[ampName].update(rowBase)
            rows[ampName]["amplifier"] = ampName
            for key, value in stats.items():
                rows[ampName][f"{self.config.stageName}_{key}"] = value

        # VERIFY results
        for ampName, stats in statisticsDict["VERIFY"]["AMP"].items():
            for key, value in stats.items():
                rows[ampName][f"{self.config.stageName}_VERIFY_{key}"] = value

        # METADATA results
        # DET results
        rows['detector'] = rowBase
        for testName, value in statisticsDict["DET"].items():
            verifyDict = statisticsDict["VERIFY"]["DET"]
            rows['detector'][f"{self.config.stageName}_DET_{testName}"] = value
            if testName in verifyDict:
                rows['detector'][f"{self.config.stageName}_DET_VERIFY_{testName}"] = verifyDict[testName]

        # ISR results
        if self.config.useIsrStatistics and "ISR" in statisticsDict:
            for ampName, stats in statisticsDict["ISR"]["CALIBDIST"].items():
                for level in self.config.expectedDistributionLevels:
                    key = f"LSST CALIB {self.config.stageName.upper()} {ampName} DISTRIBUTION {level}-PCT"
                    rows[ampName][f"{self.config.stageName}_FLAT_DIST_{level}_PCT"] = stats[key]

        # pack final list
        for ampName, stats in rows.items():
            rowList.append(stats)

        return rowList, matrixRowList


class CpVerifyFlatExpMergeConfig(CpVerifyExpMergeByFilterConfig):
    """Inherits from base CpVerifyExpMergeConfig
    """

    def setDefaults(self):
        super().setDefaults()
        self.statKeywords = {
            'EXPOSURE_SCATTER': 'STDEV',  # noqa F841
        }


class CpVerifyFlatExpMergeTask(CpVerifyExpMergeByFilterTask):
    """Inherits from base CpVerifyExpMergeTask
    """
    ConfigClass = CpVerifyFlatExpMergeConfig
    _DefaultName = 'cpVerifyFlatExpMerge'

    def calcStatistics(self, statisticsDictionary):
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

        return {'SCATTER': float(np.std(detectorMeans))}

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

    def pack(self, statisticsDict, dimensions, outKey):
        """Repack information into flat tables.

        This method should be redefined in subclasses.

        Parameters
        ----------
        statisticsDictionary : `dict` [`str`, `dict` [`str`, scalar]],
            Dictionary of measured statistics.  The inner dictionary
            should have keys that are statistic names (`str`) with
            values that are some sort of scalar (`int` or `float` are
            the mostly likely types).

        Returns
        -------
        outputResults : `list` [`dict`]
            A list of rows to add to the output table.
        outputMatrix : `list` [`dict`]
            A list of rows to add to the output matrix.
        """
        rowList = []
        matrixRowList = None

        # We can only do stats if we only have one thing.
        rowBase = {
            "instrument": dimensions[0]["instrument"],
            "exposure": dimensions[0]["exposure"],
            "detector": dimensions[0]["detector"],
        }

        # This only needs to add the new results.
        stats = statisticsDict[outKey]
        verify = statisticsDict["VERIFY"][outKey]

        for test, value in stats.items():
            rowBase[f"{self.config.stageName}_{test}"] = value
            rowBase[f"{self.config.stageName}_VERIFY_{test}"] = verify[test]

        # pack final list
        rowList.append(rowBase)

        return rowList, matrixRowList
