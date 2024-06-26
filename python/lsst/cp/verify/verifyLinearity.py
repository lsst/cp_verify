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

__all__ = [
    "CpVerifyLinearityConnections",
    "CpVerifyLinearityConfig",
    "CpVerifyLinearityTask",
    "CpvLinearitySolveConfig",
    "CpvLinearitySolveTask",
]

import lsst.cp.pipe as cpPipe
import lsst.pipe.base as pipeBase

import numpy as np
import lsst.pex.config as pexConfig

import lsst.pipe.base.connectionTypes as cT
from .verifyCalib import CpVerifyCalibConfig, CpVerifyCalibTask, CpVerifyCalibConnections


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

    maxResidualThresholdSpline = pexConfig.Field(
        dtype=float,
        doc="Maximum percentage for linearity residuals, if spline",
        default=1.0,
    )
    expectedQuadraticCoeffPolynomial = pexConfig.Field(
        dtype=float,
        doc="Expected amplitude of second-order non-linearity coefficient, if polynomial.",
        default=1e-6,
    )
    maxResidualThresholdTable = pexConfig.Field(
        dtype=float,
        doc="Maximum percentage for linearity residuals, if lookup table.",
        default=1.0,
    )

    def setDefaults(self):
        super().setDefaults()
        self.stageName = 'LINEARITY'


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
        calibMetadata = inputCalib.getMetadata()
        detId = calibMetadata['DETECTOR']
        detector = camera[detId]
        outputStatistics = {amp.getName(): {} for amp in detector}
        for amp in detector:
            ampName = amp.getName()
            outputStatistics[ampName]['FIT_PARAMS'] = inputCalib.fitParams[ampName].tolist()
            outputStatistics[ampName]['FIT_PARAMS_ERR'] = inputCalib.fitParamsErr[ampName].tolist()
            outputStatistics[ampName]['FIT_RESIDUALS'] = inputCalib.fitResiduals[ampName].tolist()
            outputStatistics[ampName]['LINEAR_FIT'] = inputCalib.linearFit[ampName].tolist()
            outputStatistics[ampName]['LINEARITY_COEFFS'] = inputCalib.linearityCoeffs[ampName].tolist()
            outputStatistics[ampName]['LINEARITY_TYPE'] = inputCalib.linearityType[ampName]
            if inputCalib.linearityType[ampName] == 'LookupTable':
                outputStatistics[ampName]['TABLE_DATA'] = inputCalib.tableData[ampName].tolist()
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
        calibMetadata = calib.getMetadata()
        detId = calibMetadata['DETECTOR']
        detector = camera[detId]
        verifyAmpStats = {}
        for amp in detector:
            verify = {}
            ampName = amp.getName()
            linearityType = calib.linearityType[ampName]
            measuredCoeffs = calib.linearityCoeffs[ampName]
            if linearityType == 'Spline':
                binCenters, values = np.split(calib.linearityCoeffs[ampName], 2)
                maxError = max(abs(values/binCenters))*100
                verify['MAX_RESIDUAL_ERROR'] = bool(maxError <= self.config.maxResidualThresholdSpline)
            elif linearityType == 'Squared':
                c0 = np.abs(measuredCoeffs[2])
                verify['MAX_RESIDUAL_ERROR'] = bool(c0 <= self.config.expectedQuadraticCoeffPolynomial)
            elif linearityType == 'Polynomial':
                epsilon = np.sqrt(self.config.expectedQuadraticCoeffPolynomial)
                coeffs = np.abs(measuredCoeffs[2:])
                # coeffs[0] is now the quadratic term. Scale higher-order terms
                # with epsilon.
                thresholds = [coeffs[0]*epsilon**order for order, c in enumerate(coeffs)]
                thresholds[0] = epsilon**2
                verify['MAX_RESIDUAL_ERROR'] = bool(np.all(coeffs <= thresholds))
            elif linearityType == 'LookupTable':
                # If 'LookupTable', linearityCoeffs is of the form {'C10':
                # array([0, 0]), 'C11': array([1, 0]), ... }
                indexTableAmp, offset = calib.linearityCoeffs[ampName]
                # Indices of second axis of table is flux range, up to
                # 2**18 ADU
                indices = np.arange(1, len(calib.tableData.T)) + offset
                # Look at correction (what the table provides) divided by
                # signal
                delta = calib.tableData[indexTableAmp, :][1:] / indices
                maxError = np.max(np.abs(delta))
                verify['MAX_RESIDUAL_ERROR'] = bool(maxError <= self.config.maxResidualThresholdTable)
            else:
                # 'None' type found. Dummy linearizer.
                self.log.warning("Dummy linearizer found (type: `None`).")
                verify['MAX_RESIDUAL_ERROR'] = False

            # Overall success among all tests for this amp.
            verify['SUCCESS'] = bool(np.all(list(verify.values())))
            if verify['SUCCESS'] is False:
                success = False

            verifyAmpStats[ampName] = verify

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
        return {'AMP': verifyAmpStats,
                'DET': verifyDetStatsFinal}, bool(success)

    def repackStats(self, statisticsDict, dimensions):
        # docstring inherited
        rows = {}
        rowList = []
        matrixRowList = []

        rowBase = {
            "instrument": dimensions["instrument"],
            "detector": dimensions["detector"],
        }

        # AMP results
        for ampName, stats in statisticsDict["AMP"].items():
            coeffs = stats.pop("LINEARITY_COEFFS")
            if stats["LINEARITY_TYPE"] == 'spline':
                centers, values = np.split(coeffs, 2)
                centers = centers.tolist()
                values = values.tolist()
            else:
                values = coeffs
                centers = np.full_like(values, np.nan).tolist()

            rows[ampName] = {}
            rows[ampName].update(rowBase)
            rows[ampName]["amplifier"] = ampName
            rows[ampName][f"{self.config.stageName}_CENTERS"] = centers
            rows[ampName][f"{self.config.stageName}_COEFFS"] = values
            for testName, value in stats.items():
                rows[ampName][f"{self.config.stageName}_{testName}"] = value

        # VERIFY results
        rows["detector"] = rowBase
        for ampName, stats in statisticsDict["VERIFY"]["AMP"].items():
            for testName, value in stats.items():
                rows[ampName][f"{self.config.stageName}_{testName}"] = value

        for ampName, stats in rows.items():
            rowList.append(stats)

        return rowList, matrixRowList


# Subclass the linearity classess so that the linearizer
# is a regular Input instead of a PrerequisiteInput
class CpvLinearitySolveConnections(pipeBase.PipelineTaskConnections,
                                   dimensions=("instrument", "detector")):
    dummy = cT.Input(
        name="raw",
        doc="Dummy exposure.",
        storageClass='Exposure',
        dimensions=("instrument", "exposure", "detector"),
        multiple=True,
        deferLoad=True,
    )
    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera Geometry definition.",
        storageClass="Camera",
        dimensions=("instrument", ),
        isCalibration=True,
    )
    inputPtc = cT.Input(
        name="ptc",
        doc="Input PTC dataset.",
        storageClass="PhotonTransferCurveDataset",
        dimensions=("instrument", "detector"),
        isCalibration=True,
    )
    inputPhotodiodeCorrection = cT.Input(
        name="cpvPdCorrection",
        doc="Input photodiode correction.",
        storageClass="IsrCalib",
        dimensions=("instrument", ),
        isCalibration=True,
    )

    outputLinearizer = cT.Output(
        name="cpvLinearity",
        doc="Output linearity measurements.",
        storageClass="Linearizer",
        dimensions=("instrument", "detector"),
        isCalibration=True,
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if config.applyPhotodiodeCorrection is not True:
            self.inputs.discard("inputPhotodiodeCorrection")


class CpvLinearitySolveConfig(cpPipe.LinearitySolveConfig,
                              pipelineConnections=CpvLinearitySolveConnections):
    pass


class CpvLinearitySolveTask(cpPipe.LinearitySolveTask):

    ConfigClass = CpvLinearitySolveConfig
    _DefaultName = "cpvLinearityTask"

    pass


class CpvPhotodiodeCorrectionConnections(pipeBase.PipelineTaskConnections,
                                         dimensions=("instrument",)):

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera Geometry definition.",
        storageClass="Camera",
        dimensions=("instrument", ),
        isCalibration=True,
    )

    inputPtc = cT.Input(
        name="ptc",
        doc="Input PTC dataset.",
        storageClass="PhotonTransferCurveDataset",
        dimensions=("instrument", "detector"),
        multiple=True,
        isCalibration=True,
    )

    inputLinearizer = cT.Input(
        name="unCorrectedLinearizer",
        doc="Raw linearizers that have not been corrected.",
        storageClass="Linearizer",
        dimensions=("instrument", "detector"),
        multiple=True,
        isCalibration=True,
    )

    outputPhotodiodeCorrection = cT.Output(
        name="pdCorrection",
        doc="Correction of photodiode systematic error.",
        storageClass="IsrCalib",
        dimensions=("instrument", ),
        isCalibration=True,
    )


class CpvPhotodiodeCorrectionConfig(cpPipe.PhotodiodeCorrectionConfig,
                                    pipelineConnections=CpvPhotodiodeCorrectionConnections):
    pass


class CpvPhotodiodeCorrectionTask(cpPipe.PhotodiodeCorrectionTask):

    ConfigClass = CpvPhotodiodeCorrectionConfig
    _DefaultName = "cpvPdCorrTask"

    pass
