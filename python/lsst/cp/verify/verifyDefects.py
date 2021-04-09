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
import scipy

from .verifyStats import CpVerifyStatsConfig, CpVerifyStatsTask, CpVerifyStatsConnections


__all__ = ['CpVerifyDefectsConfig', 'CpVerifyDefectsTask']


class CpVerifyDefectsConfig(CpVerifyStatsConfig,
                            pipelineConnections=CpVerifyStatsConnections):
    """Inherits from base CpVerifyStatsConfig
    """

    def setDefaults(self):
        super().setDefaults()
        self.maskNameList = ['BAD']  # noqa F821

        self.imageStatKeywords = {'DEFECT_PIXELS': 'NMASKED',  # noqa F821
                                  'OUTLIERS': 'NCLIPPED',
                                  'MEDIAN': 'MEDIAN',
                                  'STDEV': 'STDEVCLIP',
                                  'MIN': 'MIN',
                                  'MAX': 'MAX', }
        self.unmaskedImageStatKeywords = {'UNMASKED_MIN': 'MIN',  # noqa F821
                                          'UNMASKED_MAX': 'MAX',
                                          'UNMASKED_OUTLIERS': 'NCLIPPED', }


class CpVerifyDefectsTask(CpVerifyStatsTask):
    """Defects verification sub-class, implementing the verify method.

    This also applies additional image processing statistics.
    """
    ConfigClass = CpVerifyDefectsConfig
    _DefaultName = 'cpVerifyDefects'

    def imageStatistics(self, exposure, statControl):
        """Measure additional defect statistics.

        This calls the parent class method first, then adds additional
        measurements.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure containing the ISR processed data to measure.
        statControl : `lsst.afw.math.StatControl`
            Statistics control object with the parameters defined by
            the config.

        Returns
        -------
        outputStats : `dict` [`str`, `dict` [`str`, scalar]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the statistics measured and their values.
        """
        outputStats = super().imageStatistics(exposure, statControl)

        # Is this a useful test?  It saves having to do chi^2 fits,
        # which are going to be biased by the bulk of points.
        for amp in exposure.getDetector():
            ampName = amp.getName()
            ampExp = exposure.Factory(exposure, amp.getBBox())

            normImage = ampExp.getImage()
            normArray = normImage.getArray()

            normArray -= outputStats[ampName]['MEDIAN']
            normArray /= outputStats[ampName]['STDEV']

            probability = scipy.stats.norm.pdf(normArray)
            trials = np.prod(probability.shape)
            outliers = np.where(probability < 1.0 / trials, 1.0, 0.0)
            outputStats[ampName]['STAT_OUTLIERS'] = np.sum(outliers)

        return outputStats

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
        ampStats = statisticsDictionary['AMP']
        verifyStats = dict()
        success = True
        for ampName, stats in ampStats.items():
            verify = dict()

            # These are not defined yet.
            if stats['UNMASKED_OUTLIERS'] > stats['OUTLIERS']:
                verify['OUTLIERS'] = True
            else:
                verify['OUTLIERS'] = False

            verify['MIN'] = True if stats['MIN'] >= stats['UNMASKED_MIN'] else False
            verify['MAX'] = True if stats['MAX'] <= stats['UNMASKED_MAX'] else False

            verify['PROB_TEST'] = True if stats['STAT_OUTLIERS'] == stats['DEFECT_PIXELS'] else False

            verify['SUCCESS'] = bool(np.all(list(verify.values())))
            if verify['SUCCESS'] is False:
                success = False

            # Why is this necessary?  Why are these booleans not
            # boolean enough for yaml.safe_dump?
            verify = {key: bool(value) for key, value in verify.items()}

            verifyStats[ampName] = verify

        return {'AMP': verifyStats}, bool(success)
