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
from scipy.optimize import least_squares

from .verifyCalib import CpVerifyCalibConfig, CpVerifyCalibTask, CpVerifyCalibConnections

__all__ = ['CpVerifyPtcConfig', 'CpVerifyPtcTask']


class CpVerifyPtcConfig(CpVerifyCalibConfig,
                        pipelineConnections=CpVerifyCalibConnections):
    """Inherits from base CpVerifyCalibConfig."""

    def setDefaults(self):
        super().setDefaults()

    gainThreshold = pexConfig.Field(
        dtype=float,
        doc="Maximum percentage difference between PTC gain and nominal amplifier gain.",
        default=5.0,
    )

    noiseThreshold = pexConfig.Field(
        dtype=float,
        doc="Maximum percentage difference between PTC readout noise and nominal "
            "amplifier readout noise.",
        default=5.0,
    )

    turnoffThreshold = pexConfig.Field(
        dtype=float,
        doc="Minimun full well requirement (in electrons). To be compared with the "
            "reported PTC turnoff per amplifier.",
        default=90000,
    )

    a00MinITL = pexConfig.Field(
        dtype=float,
        doc="Minimum a00 (c.f., Astier+19) for ITL CCDs.",
        default=-4.56e-6,
    )

    a00MaxITL = pexConfig.Field(
        dtype=float,
        doc="Maximum a00 (c.f., Astier+19) for ITL CCDs.",
        default=6.91e-7,
    )

    a00MinE2V = pexConfig.Field(
        dtype=float,
        doc="Minimum a00 (c.f., Astier+19) for E2V CCDs.",
        default=-3.52e-6,
    )

    a00MaxE2V = pexConfig.Field(
        dtype=float,
        doc="Maximum a00 (c.f., Astier+19) for E2V CCDs.",
        default=-2.61e-6,
    )


def linearModel(x, m, b):
    """A linear model.
    """
    return m*x + b


def modelResidual(p, x, y):
    """Model residual for fit below.
    """
    return y - linearModel(x, *p)


class CpVerifyPtcTask(CpVerifyCalibTask):
    """PTC verification sub-class, implementing the verify method.
    """
    ConfigClass = CpVerifyPtcConfig
    _DefaultName = 'cpVerifyPtc'

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
        ptcFitType = calibMetadata['PTC_FIT_TYPE']
        outputStatistics = {amp.getName(): {} for amp in detector}
        for amp in detector:
            ampName = amp.getName()
            calibGain = inputCalib.gain[ampName]
            outputStatistics[ampName]['PTC_GAIN'] = calibGain
            outputStatistics[ampName]['AMP_GAIN'] = amp.getGain()
            outputStatistics[ampName]['PTC_NOISE'] = inputCalib.noise[ampName]
            outputStatistics[ampName]['AMP_NOISE'] = amp.getReadNoise()
            outputStatistics[ampName]['PTC_TURNOFF'] = inputCalib.ptcTurnoff[ampName]
            outputStatistics[ampName]['PTC_FIT_TYPE'] = ptcFitType
            outputStatistics[ampName]['PTC_ROW_MEAN_VARIANCE'] = inputCalib.rowMeanVariance[ampName].tolist()
            outputStatistics[ampName]['PTC_MAX_RAW_MEANS'] = float(np.max(inputCalib.rawMeans[ampName]))
            # To plot Covs[ij] vs flux
            rawFlux = inputCalib.rawMeans[ampName].tolist()
            outputStatistics[ampName]['PTC_RAW_MEANS'] = rawFlux
            mask = inputCalib.expIdMask[ampName].tolist()
            outputStatistics[ampName]['PTC_EXP_ID_MASK'] = mask
            covs = inputCalib.covariances[ampName]
            outputStatistics[ampName]['PTC_COV_10'] = covs[:, 1, 0].tolist()
            outputStatistics[ampName]['PTC_COV_01'] = covs[:, 0, 1].tolist()
            outputStatistics[ampName]['PTC_COV_11'] = covs[:, 1, 1].tolist()
            outputStatistics[ampName]['PTC_COV_20'] = covs[:, 2, 0].tolist()
            outputStatistics[ampName]['PTC_COV_02'] = covs[:, 0, 2].tolist()
            # Calculate and save the slopes and offsets from Covs[ij] vs flux
            keys = ['PTC_COV_10', 'PTC_COV_01', 'PTC_COV_11', 'PTC_COV_20',
                    'PTC_COV_02']
            maskedFlux = np.array(rawFlux)[mask]
            for key in keys:
                maskedCov = np.array(outputStatistics[ampName][key])[mask]
                linearFit = least_squares(modelResidual, [1., 0.0],
                                          args=(np.array(maskedFlux), np.array(maskedCov)),
                                          loss='cauchy')
                slopeKey = key + '_FIT_SLOPE'
                offsetKey = key + '_FIT_OFFSET'
                successKey = key + '_FIT_SUCCESS'
                outputStatistics[ampName][slopeKey] = float(linearFit.x[0])
                outputStatistics[ampName][offsetKey] = float(linearFit.x[1])
                outputStatistics[ampName][successKey] = linearFit.success

            if ptcFitType == 'EXPAPPROXIMATION':
                outputStatistics[ampName]['PTC_BFE_A00'] = float(inputCalib.ptcFitPars[ampName][0])
            if ptcFitType == 'FULLCOVARIANCE':
                outputStatistics[ampName]['PTC_BFE_A00'] = float(inputCalib.aMatrix[ampName][0][0])

            # Test from eo_pipe: github.com/lsst-camera-dh/eo-pipe;
            # ptcPlotTask.py
            # Slope of [variance of means of rows](electrons^2)
            # vs [2*signal(electrons)/numCols]
            numCols = amp.getBBox().width
            mask = inputCalib.expIdMask[ampName]
            rowMeanVar = inputCalib.rowMeanVariance[ampName][mask]*calibGain**2
            signal = inputCalib.rawMeans[ampName][mask]*calibGain
            try:
                slope = sum(rowMeanVar) / sum(2.*signal/numCols)
            except ZeroDivisionError:
                slope = np.nan
            outputStatistics[ampName]['ROW_MEAN_VARIANCE_SLOPE'] = float(slope)

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
        ptcFitType = calibMetadata['PTC_FIT_TYPE']
        # 'DET_SER' is of the form 'ITL-3800C-229'
        detVendor = calibMetadata['DET_SER'].split('-')[0]

        for amp in detector:
            verify = {}
            ampName = amp.getName()
            calibGain = calib.gain[ampName]

            diffGain = (np.abs(calibGain - amp.getGain()) / amp.getGain())*100
            diffNoise = (np.abs(calib.noise[ampName] - amp.getReadNoise()) / amp.getReadNoise())*100

            # DMTN-101: 16.1 and 16.2
            # The fractional relative difference between the fitted PTC and the
            # nominal amplifier gain and readout noise values should be less
            # than a certain threshold (default: 5%).
            verify['PTC_GAIN'] = bool(diffGain < self.config.gainThreshold)
            verify['PTC_NOISE'] = bool(diffNoise < self.config.noiseThreshold)

            # DMTN-101: 16.3
            # Check that the measured PTC turnoff is at least greater than the
            # full-well requirement of 90k e-.
            turnoffCut = self.config.turnoffThreshold
            verify['PTC_TURNOFF'] = bool(calib.ptcTurnoff[ampName]*calibGain > turnoffCut)
            # DMTN-101: 16.4
            # Check the a00 value (brighter-fatter effect).
            # This is a purely electrostatic parameter that should not change
            # unless voltages are changed (e.g., parallel, bias voltages).
            # Check that the fitted a00 parameter per CCD vendor is within a
            # range motivated by measurements on data (DM-30171).
            if ptcFitType in ['EXPAPPROXIMATION', 'FULLCOVARIANCE']:
                # a00 is a fit parameter from these models.
                if ptcFitType == 'EXPAPPROXIMATION':
                    a00 = calib.ptcFitPars[ampName][0]
                else:
                    a00 = calib.aMatrix[ampName][0][0]
                if detVendor == 'ITL':
                    a00Max = self.config.a00MaxITL
                    a00Min = self.config.a00MinITL
                    verify['PTC_BFE_A00'] = bool(a00 > a00Min and a00 < a00Max)
                elif detVendor == 'E2V':
                    a00Max = self.config.a00MaxE2V
                    a00Min = self.config.a00MinE2V
                    verify['PTC_BFE_A00'] = bool(a00 > a00Min and a00 < a00Max)
                else:
                    raise RuntimeError(f"Detector type {detVendor} not one of 'ITL' or 'E2V'")

            # Overall success among all tests for this amp.
            verify['SUCCESS'] = bool(np.all(list(verify.values())))
            if verify['SUCCESS'] is False:
                success = False

            verifyStats[ampName] = verify

        return {'AMP': verifyStats}, bool(success)
