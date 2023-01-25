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
import lsst.pex.config as pexConfig

import lsst.pipe.base.connectionTypes as cT
from .verifyCalib import CpVerifyCalibConfig, CpVerifyCalibTask, CpVerifyCalibConnections

__all__ = ['CpVerifyLinearityConnections', 'CpVerifyLinearityConfig', 'CpVerifyLinearityTask']


class CpVerifyLinearityConnections(CpVerifyCalibConnections,
                                   dimensions={"instrument", "detector"},
                                   defaultTemplates={}):
    inputCalib = cT.Input(
        name="calib",
        doc="Input calib to calculate statistics for.",
        storageClass="Linearizer",
        dimensions=("instrument", "detector"),
        isCalibration=True
    )


class CpVerifyLinearityConfig(CpVerifyCalibConfig,
                              pipelineConnections=CpVerifyLinearityConnections):
    """Inherits from base CpVerifyCalibConfig."""

    def setDefaults(self):
        super().setDefaults()

    maxResidualThreshold = pexConfig.Field(
        dtype=float,
        doc="Maximum percentage for linearity residuals.",
        default=1.0,
    )


class CpVerifyLinearityTask(CpVerifyCalibTask):
    """Linearity verification sub-class, implementing the verify method.
    """
    ConfigClass = CpVerifyLinearityConfig
    _DefaultName = 'cpVerifyLinearity'

    def detectorStatistics(self, inputCalib, camera=None):
        """Calculate detector level statistics from the calibration.

        Parameters
        ----------
        inputCalib : `lsst.ip.isr.IsrCalib`
            The calibration to verify.
        camera : `lsst.afw.cameraGeom.Camera`, optional
            Input camera to get detectors from.

        Returns
        -------
        outputStatistics : `dict` [`str`, scalar]
            A dictionary of the statistics measured and their values.
        """
        return {}

    def amplifierStatistics(self, inputCalib, camera=None):
        """Calculate detector level statistics from the calibration.

        Parameters
        ----------
        inputCalib : `lsst.ip.isr.IsrCalib`
            The calibration to verify.
        camera : `lsst.afw.cameraGeom.Camera`, optional
            Input camera to get detectors from.

        Returns
        -------
        outputStatistics : `dict` [`str`, scalar]
            A dictionary of the statistics measured and their values.
        """
        calibMetadata = inputCalib.getMetadata().toDict()
        detId = calibMetadata['DETECTOR']
        detector = camera[detId]
        outputStatistics = {amp.getName(): {} for amp in detector}
        for amp in detector:
            ampName = amp.getName()
            outputStatistics[ampName]['FIT_CHI_SQ'] = inputCalib.fitChiSq
            outputStatistics[ampName]['FIT_PARAMS'] = inputCalib.fitParams
            outputStatistics[ampName]['FIT_PARAMS_ERR'] = inputCalib.fitParamsErr
            outputStatistics[ampName]['FIT_RESIDUALS'] = inputCalib.fitResiduals
            outputStatistics[ampName]['LINEAR_FIT'] = inputCalib.linearFit[ampName]
            outputStatistics[ampName]['LINEARITY_COEFFS'] = inputCalib.linearityCoeffs[ampName]
            outputStatistics[ampName]['LINEARITY_TYPE'] = inputCalib.linearityType[ampName]
            outputStatistics[ampName]['TABLE_DATA'] = inputCalib.tableData[ampName]

        return outputStatistics

    def verify(self, calib, statisticsDict, camera=None):
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
        camera : `lsst.afw.cameraGeom.Camera`, optional
            Input camera to get detectors from.

        Returns
        -------
        outputStatistics : `dict` [`str`, `dict` [`str`, `bool`]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the verification criteria.
        success : `bool`
            A boolean indicating whether all tests have passed.
        """
        verifyStats = {}
        success = True
        calibMetadata = calib.getMetadata().toDict()
        detId = calibMetadata['DETECTOR']
        detector = camera[detId]

        for amp in detector:
            verify = {}
            ampName = amp.getName()
            linearityType = calib.linearityType[ampName]
            if linearityType == 'Spline':
                binCenters, values = np.split(calib.linearityCoeffs[ampName], 2)
                maxError = maxError = max(abs(values/binCenters))*100
            elif linearityType == 'Squared':
                # Define something for the other types
                maxError = 100
            elif linearityType == 'Polynomial':
                maxError = 100
            elif linearityType == 'LookupTable':
                maxError = 100
            else:
                # 'None' "Create a dummy solution.
                maxError = 100

            verify['MAX_RESIDUAL_ERROR'] = bool(maxError < self.config.maxResidualThreshold)

            # Overall success among all tests for this amp.
            verify['SUCCESS'] = bool(np.all(list(verify.values())))
            if verify['SUCCESS'] is False:
                success = False

            verifyStats[ampName] = verify

        # Loop over amps to make a detector summary.
        verifyDetStats = {'MAX_RESIDUAL_ERROR': []}
        for amp in verifyStats:
            for testName in verifyStats[amp]:
                if testName == 'SUCCESS':
                    continue
                verifyDetStats[testName].append(verifyStats[amp][testName])

        # VerifyDetStatsFinal has final boolean test over all amps
        verifyDetStatsFinal = {}
        for testName in verifyDetStats:
            testBool = bool(np.all(list(verifyDetStats[testName])))
            # Save the tests that failed
            if not testBool:
                verifyDetStatsFinal[testName] = bool(np.all(list(verifyDetStats[testName])))

        return verifyDetStatsFinal, bool(success)
