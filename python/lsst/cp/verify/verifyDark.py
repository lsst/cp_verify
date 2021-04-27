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
        self.imageStatKeywords = {'MEAN': 'MEAN',  # noqa F841
                                  'NOISE': 'STDEVCLIP', }
        self.crImageStatKeywords = {'CR_NOISE': 'STDEV', }  # noqa F841


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
        detector = exposure.getDetector()
        ampStats = statisticsDict['AMP']
        verifyStats = {}
        success = True
        for ampName, stats in ampStats.items():
            verify = {}

            # DMTN-101 Test 5.2: Mean is 0.0:
            verify['MEAN'] = bool(np.abs(stats['MEAN']) < stats['NOISE'])

            # DMTN-101 Test 5.3: Clipped mean matches readNoise
            amp = detector[ampName]
            verify['NOISE'] = bool(np.abs(stats['NOISE'] - amp.getReadNoise())/amp.getReadNoise() <= 0.05)

            # DMTN-101 Test 5.4: CR rejection matches clipped mean
            verify['CR_NOISE'] = bool(np.abs(stats['NOISE'] - stats['CR_NOISE'])/stats['CR_NOISE'] <= 0.05)

            verify['SUCCESS'] = bool(np.all(list(verify.values())))
            if verify['SUCCESS'] is False:
                success = False

            verifyStats[ampName] = verify

        return {'AMP': verifyStats}, bool(success)
