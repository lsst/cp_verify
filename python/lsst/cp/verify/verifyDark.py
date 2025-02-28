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

from .verifyStats import CpVerifyStatsConfig, CpVerifyStatsTask, CpVerifyStatsConnections

__all__ = ['CpVerifyDarkConfig', 'CpVerifyDarkTask']


class CpVerifyDarkConfig(CpVerifyStatsConfig,
                         pipelineConnections=CpVerifyStatsConnections):
    """Inherits from base CpVerifyStatsConfig.
    """

    def setDefaults(self):
        super().setDefaults()
        self.stageName = 'DARK'
        self.imageStatKeywords = {'MEAN': 'MEAN',  # noqa F841
                                  'NOISE': 'STDEVCLIP', }
        self.crImageStatKeywords = {'CR_NOISE': 'STDEV', }  # noqa F841
        self.metadataStatKeywords = {
            'LSST ISR OVERSCAN RESIDUAL SERIAL STDEV': 'READ_NOISE',
        }  # noqa F841


class CpVerifyDarkTask(CpVerifyStatsTask):
    """Dark verification sub-class, implementing the verify method.
    """
    ConfigClass = CpVerifyDarkConfig
    _DefaultName = 'cpVerifyDark'

    def verify(self, exposure, statisticsDict):
        """Verify that the measured statistics meet the verification criteria.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure the statistics are from.
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
        ampStats = statisticsDict['AMP']
        metadataStats = statisticsDict['METADATA']

        verifyStats = {}
        success = True
        for ampName, stats in ampStats.items():
            verify = {}

            # DMTN-101 Test 5.2: Mean is 0.0:
            verify['MEAN'] = bool(np.abs(stats['MEAN']) < stats['NOISE'])

            # DMTN-101 Test 5.3: Clipped mean matches readNoise.  This
            # test should use the nominal detector read noise.  The
            # f"READ_NOISE" metadata entry contains the
            # measured dispersion in the overscan-corrected overscan
            # region, which should provide an estimate of the read
            # noise.  However, directly using this value will cause
            # some fraction of verification runs to fail if the
            # scatter in read noise values is comparable to the test
            # threshold, as the overscan residual measured may be
            # sampling from the low end tail of the distribution.
            # This measurement is also likely to be smaller than that
            # measured on the bulk of the image as the overscan
            # correction should be an optimal fit to the overscan
            # region, but not necessarily for the image region.
            readNoise = exposure.getMetadata().get(f"LSST ISR READNOISE {ampName}", np.nan)
            verify['NOISE'] = bool((stats['NOISE'] - readNoise)/readNoise <= 0.05)

            # DMTN-101 Test 5.4: CR rejection matches clipped mean
            verify['CR_NOISE'] = bool(np.abs(stats['NOISE'] - stats['CR_NOISE'])/stats['CR_NOISE'] <= 0.05)

            # Confirm this hasn't triggered a raise condition.
            if 'FORCE_FAILURE' in stats:
                verify['PROCESSING'] = False

            verify['SUCCESS'] = bool(np.all(list(verify.values())))
            if verify['SUCCESS'] is False:
                success = False

            # After determining the verification status for this
            # exposure, we can also check to see how well the read
            # noise measured from the overscan residual matches the
            # nominal value used above in Test 5.3.  If these disagree
            # consistently and significantly, then the assumptions
            # used in that test may be incorrect, and the nominal read
            # noise may need recalculation.  Only perform this check
            # if the metadataStats contain the required entry.  This
            # is in ADU (the serial overscan is measured prior to gain
            # normalization), so we need to convert to electrons here.
            gain = exposure.getMetadata().get(f"LSST ISR GAIN {ampName}", np.nan)
            overscanReadNoise = gain * metadataStats['READ_NOISE'][ampName]
            if overscanReadNoise:
                if ((overscanReadNoise - readNoise)/readNoise > 0.05) or not np.isfinite(overscanReadNoise):
                    verify['READ_NOISE_CONSISTENT'] = False

            verifyStats[ampName] = verify

        return {'AMP': verifyStats}, bool(success)

    def repackStats(self, statisticsDict, dimensions):
        # docstring inherited
        rows = {}
        rowList = []
        matrixRowList = None

        if self.config.useIsrStatistics:
            mjd = statisticsDict["ISR"]["MJD"]
        else:
            mjd = np.nan

        rowBase = {
            "instrument": dimensions["instrument"],
            "exposure": dimensions["exposure"],
            "detector": dimensions["detector"],
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
        for ampName, value in statisticsDict["METADATA"]["READ_NOISE"].items():
            rows[ampName][f"{self.config.stageName}_READ_NOISE"] = value

        # ISR results
        if self.config.useIsrStatistics and "ISR" in statisticsDict:
            for ampName, stats in statisticsDict["ISR"]["CALIBDIST"].items():
                for level in self.config.expectedDistributionLevels:
                    key = f"LSST CALIB {self.config.stageName.upper()} {ampName} DISTRIBUTION {level}-PCT"
                    rows[ampName][f"{self.config.stageName}_DARK_DIST_{level}_PCT"] = stats[key]

        # pack final list
        for ampName, stats in rows.items():
            rowList.append(stats)

        return rowList, matrixRowList
