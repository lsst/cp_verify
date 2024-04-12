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
        self.stageName = 'CROSSTALK'
        self.calibStatKeywords = {'SNR': ''}  # noqa F841


class CpVerifyCrosstalkTask(CpVerifyCalibTask):
    """Crosstalk verification sub-class, implementing the verify method.
    """
    ConfigClass = CpVerifyCrosstalkConfig
    _DefaultName = 'cpVerifyCrosstalk'

    def detectorStatistics(self, inputCalib, camera=None):
        """Calculate detector level statistics from the calibration.

        Parameters
        ----------
        inputCalib : `lsst.ip.isr.IsrCalib`
            The calibration to verify.
        camera : `lsst.afw.cameraGeom.Camera`, optional
            Input camera.

        Returns
        -------
        outputStatistics : `dict` [`str`, scalar]
            A dictionary of the statistics measured and their values.
        """
        outputStatistics = {}
        outputStatistics['N_VALID'] = int(np.sum(inputCalib.coeffValid))
        outputStatistics['N_AMP'] = inputCalib.nAmp
        # I think this is the residual set, which isn't what we want,
        # but will serve as a placeholder.
        outputStatistics['COEFFS'] = inputCalib.coeffs.tolist()

        return outputStatistics

    def amplifierStatistics(self, inputCalib, camera=None):
        """Calculate amplifier level statistics from the calibration.

        Parameters
        ----------
        inputCalib : `lsst.ip.isr.IsrCalib`
            The calibration to verify.
        camera : `lsst.afw.cameraGeom.Camera`, optional
            Input camera.

        Returns
        -------
        outputStatistics : `dict` [`str`, scalar]
            A dictionary of the statistics measured and their values.
        """

        pass

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
            Input camera.

        Returns
        -------
        outputStatistics : `dict` [`str`, `dict` [`str`, `bool`]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the verification criteria.
        success : `bool`
            A boolean indicating whether all tests have passed.
        """
        verifyStats = {}
        detectorStats = statisticsDict['DET']
        success = True
        verifyStats['NO_SIGNIFICANT_DETECTION'] = True

        if detectorStats['N_VALID'] > 0:
            verifyStats['NO_SIGNIFICANT_DETECTION'] = False
            success = False

        return verifyStats, success

    def repackStats(self, statisticsDict, dimensions):
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
        matrixRowList = []

        if self.config.useIsrStatistics:
            mjd = statisticsDict["ISR"]["MJD"]
        else:
            mjd = np.nan

        rowBase = {
            "instrument": dimensions["instrument"],
            "detector": dimensions["detector"],
            "mjd": mjd,
            "amplifier": "detector",
        }
        row = {}
        row.update(rowBase)

        # Pack DET results
        for key, value in statisticsDict['DET'].items():
            if key == 'COEFFS':
                matrixRowBase = {
                    "instrument": dimensions["instrument"],
                    "detector": dimensions["detector"],
                    "detectorComp": dimensions["detector"],
                    "mjd": mjd,
                }

                Umax = statisticsDict['DET']['N_AMP']
                Vmax = statisticsDict['DET']['N_AMP']
                for u in range(Umax):
                    for v in range(Vmax):
                        matrixRow = matrixRowBase
                        matrixRow["amplifierIdx"] = u
                        matrixRow["amplifierCompIdx"] = v
                        matrixRow["coefficient"] = value[u][v]

                        matrixRowList.append(matrixRow)
            else:
                row[f"{self.config.stageName}_{key}"] = value

        # VERIFY results

        for key, value in statisticsDict["VERIFY"].items():
            row[f"{self.config.stageName}_VERIFY_{key}"] = value

        # pack final list
        rowList.append(row)

        return rowList, matrixRowList
