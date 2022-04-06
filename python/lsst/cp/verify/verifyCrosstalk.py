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

import lsst.pipe.base.connectionTypes as cT
from .verifyCalib import CpVerifyCalibConfig, CpVerifyCalibTask, CpVerifyCalibConnections


__all__ = ['CpVerifyCrosstalkConfig', 'CpVerifyCrosstalkTask']


class CpVerifyCrosstalkConnections(CpVerifyCalibConnections,
                                   dimensions={"instrument", "detector"},
                                   defaultTemplates={}):
    inputCalib = cT.Input(
        name="calib",
        doc="Input calib to calculate statistics for.",
        storageClass="CrosstalkCalib",
        dimensions=["instrument", "detector"],
        isCalibration=True
    )


class CpVerifyCrosstalkConfig(CpVerifyCalibConfig,
                              pipelineConnections=CpVerifyCrosstalkConnections):
    """Inherits from base CpVerifyCalibConfig."""

    def setDefaults(self):
        super().setDefaults()
        self.calibStatKeywords = {'SNR': ''}  # noqa F841


class CpVerifyCrosstalkTask(CpVerifyCalibTask):
    """Crosstalk verification sub-class, implementing the verify method.
    """
    ConfigClass = CpVerifyCrosstalkConfig
    _DefaultName = 'cpVerifyCrosstalk'

    def detectorStatistics(self, inputCalib):
        """Calculate detector level statistics from the calibration.

        Parameters
        ----------
        inputCalib : `lsst.ip.isr.IsrCalib`
            The calibration to verify.

        Returns
        -------
        outputStatistics : `dict` [`str`, scalar]
            A dictionary of the statistics measured and their values.
        """
        outputStatistics = {}
        outputStatistics['N_VALID'] = int(np.sum(inputCalib.coeffValid))
        outputStatistics['N_AMP'] = inputCalib.nAmp

        return outputStatistics

    def verify(self, calib, statisticsDict):
        """Verify that the calibration meets the verification criteria.

        Parameters
        ----------
        inputCalib : `lsst.ip.isr.IsrCalib`
            The calibration to verify.
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
            A boolean indicating whether all tests have passed.
        """
        verifyStats = {}
        detectorStats = statisticsDict['DETECTOR']
        success = True
        verifyStats['NO_SIGNIFICANT_DETECTION'] = True

        if detectorStats['N_VALID'] > 0:
            verifyStats['NO_SIGNIFICANT_DETECTION'] = False
            success = False

        return verifyStats, success
