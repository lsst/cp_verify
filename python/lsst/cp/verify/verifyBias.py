#
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
#
import numpy as np

from .verifyStats import CpVerifyStatsConfig, CpVerifyStatsTask, CpVerifyStatsConnections


__all__ = ['CpVerifyBiasConfig', 'CpVerifyBiasTask']


class CpVerifyBiasConfig(CpVerifyStatsConfig,
                         pipelineConnections=CpVerifyStatsConnections):
    """Inherits from base CpVerifyStatsConfig
    """

    def setDefaults(self):
        super().setDefaults()
        self.imageStatKeywords = {'MEAN': 'MEAN',  # noqa F841
                                  'NOISE': 'STDEVCLIP', }

        self.crImageStatKeywords = {'CR_NOISE': 'STDEV', }  # noqa F841


class CpVerifyBiasTask(CpVerifyStatsTask):
    """Bias verification sub-class, implementing the verify method.
    """
    ConfigClass = CpVerifyBiasConfig
    _DefaultName = 'cpVerifyBias'

    def verify(self, exposure, statisticsDictionary):
        """Verify if the measured statistics meet the verification criteria.

        Parameters
        ----------
        statisticsDictionary : `dict` [`str`, `dict` [`str`, scalar]],
            Dictionary of measured statistics.

        Returns
        -------
        outputStats : `dict` [`str`, `dict` [`str`, `bool`]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the verification criteria.
        success : `bool`
            A boolean indicating if all tests have passed.
        """

        detector = exposure.getDetector()
        ampStats = statisticsDictionary['AMP']

        verifyStats = dict()
        success = True
        for ampName, stats in ampStats.items():
            verify = dict()

            # Test 4.2: Mean is 0.0:
            verify['MEAN'] = (np.abs(stats['MEAN']) < stats['NOISE'])

            # Test 4.3: Clipped mean matches readNoise
            amp = detector[ampName]
            verify['NOISE'] = (np.abs(stats['NOISE'] - amp.getReadNoise())/amp.getReadNoise() <= 0.05)

            # Test 4.4: CR rejection matches clipped mean
            verify['CR_NOISE'] = (np.abs(stats['NOISE'] - stats['CR_NOISE'])/stats['CR_NOISE'] <= 0.05)

            verify['SUCCESS'] = bool(np.all(list(verify.values())))
            if verify['SUCCESS'] is False:
                success = False

            # Why is this necessary?  Why are these booleans not
            # boolean enough for yaml.safe_dump?
            verify = {key: bool(value) for key, value in verify.items()}

            verifyStats[ampName] = verify

        return {'AMP': verifyStats}, bool(success)
