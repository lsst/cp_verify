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

__all__ = ['CpVerifyGainConnections', 'CpVerifyGainConfig', 'CpVerifyGainTask']


class CpVerifyGainConnections(CpVerifyCalibConnections,
                              dimensions={"instrument", "detector"},
                              defaultTemplates={}):
    exposure = cT.Input(
        name="raw",
        doc="Exposure ID of first flat from flat pair.",
        storageClass='Exposure',
        dimensions=("instrument", "detector", "exposure"),
        multiple=True,
        deferLoad=True,
    )

    inputCalib = cT.Input(
        name="calib",
        doc="Input calib to calculate statistics for.",
        storageClass="PhotonTransferCurveDataset",
        dimensions=("instrument", "detector", "exposure"),
        isCalibration=True
    )


class CpVerifyGainConfig(CpVerifyCalibConfig,
                         pipelineConnections=CpVerifyGainConnections):
    """Inherits from base CpVerifyCalibConfig."""

    gainThreshold = pexConfig.Field(
        dtype=float,
        doc="Maximum percentage difference between gain from flat pairs and nominal amplifier gain.",
        default=5.0,
    )

    noiseThreshold = pexConfig.Field(
        dtype=float,
        doc="Maximum percentage difference between overscan readout noise and nominal "
            "amplifier readout noise.",
        default=5.0,
    )

    def setDefaults(self):
        super().setDefaults()
        self.stageName = 'GAIN'


class CpVerifyGainTask(CpVerifyCalibTask):
    """Gain from flat pairs verification sub-class. Implements verify method.
    """
    ConfigClass = CpVerifyGainConfig
    _DefaultName = 'cpVerifyGain'

    def detectorStatistics(self, inputCalib, camera=None, exposure=None):
        """Calculate detector level statistics from the calibration.

        Parameters
        ----------
        inputCalib : `lsst.ip.isr.IsrCalib`
            The calibration to verify.
        camera : `lsst.afw.cameraGeom.Camera`, optional
            Input camera to get detectors from.
        exposure : `lsst.afw.image.exposure.ExposureF`, optional
            First flat-field image from pair of flats used to
            estimate the gain.

        Returns
        -------
        outputStatistics : `dict` [`str`, scalar]
            A dictionary of the statistics measured and their values.
        """
        return {}

    def amplifierStatistics(self, inputCalib, camera=None, exposure=None):
        """Calculate detector level statistics from the calibration.

        Parameters
        ----------
        inputCalib : `lsst.ip.isr.IsrCalib`
            The calibration to verify.
        camera : `lsst.afw.cameraGeom.Camera`, optional
            Input camera to get detectors from.
        exposure : `lsst.afw.image.exposure.ExposureF`, optional
            First flat-field image from pair of flats used to
            estimate the gain.

        Returns
        -------
        outputStatistics : `dict` [`str`, scalar]
            A dictionary of the statistics measured and their values.
        """
        calibMetadata = inputCalib.getMetadata().toDict()
        detId = calibMetadata['DETECTOR']
        detector = camera[detId]
        # 'DET_SER' is of the form 'ITL-3800C-229'
        detVendor = calibMetadata['DET_SER'].split('-')[0]
        # Adjust gain estimated from flat pair for flux bias, see DM-35790
        if detVendor == 'ITL':
            slope = 0.00027  # %/ADU
        elif detVendor == 'e2V':
            slope = 0.00046  # %/ADU
        outputStatistics = {amp.getName(): {} for amp in detector}

        for amp in detector:
            ampName = amp.getName()
            # Flux correction to gain-per-flat-pair method, see DM-35790
            correction = 1. - slope*inputCalib.rawMeans[ampName][0]/100
            outputStatistics[ampName]['MEAN_FLUX_FLAT_PAIR'] = inputCalib.rawMeans[ampName][0]
            outputStatistics[ampName]['GAIN_FROM_FLAT_PAIR'] = inputCalib.gain[ampName]
            outputStatistics[ampName]['GAIN_CORRECTION_FACTOR'] = correction
            outputStatistics[ampName]['AMP_GAIN'] = amp.getGain()
            outputStatistics[ampName]['ISR_NOISE'] = inputCalib.noise[ampName]
            outputStatistics[ampName]['AMP_NOISE'] = amp.getReadNoise()

        return outputStatistics

    def verify(self, calib, statisticsDict, camera=None, exposure=None):
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
        exposure : `lsst.afw.image.exposure.ExposureF`, optional
            First flat-field image from pair of flats used to
            estimate the gain.

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
        # 'DET_SER' is of the form 'ITL-3800C-229'
        detVendor = calibMetadata['DET_SER'].split('-')[0]
        # Adjust gain estimated from flat pair for flux bias, see DM-35790
        if detVendor == 'ITL':
            slope = 0.00027  # %/ADU
        elif detVendor == 'e2V':
            slope = 0.00046  # %/ADU

        for amp in detector:
            verify = {}
            ampName = amp.getName()

            # Gain from a pair of flats and noise from overscan after ISR.
            # See DM-35790.
            correction = 1. - slope*np.array(calib.rawMeans[ampName][0])/100
            gain = correction*calib.gain[ampName]

            diffGain = (np.abs(gain - amp.getGain()) / amp.getGain())*100
            diffNoise = (np.abs(calib.noise[ampName] - amp.getReadNoise()) / amp.getReadNoise())*100

            # DMTN-101: 6.1 and 6.2
            # Estimate gain from a pair of flats and compare it with the value
            # in the amplifiers.
            verify['GAIN_FROM_FLAT_PAIR'] = bool(diffGain < self.config.gainThreshold)
            # Check the empirical noise (as oppossed to fitted noise
            # from the PTC) calculated from the overscan after ISR.
            verify['ISR_NOISE'] = bool(diffNoise < self.config.noiseThreshold)

            # Overall success among all tests for this amp.
            verify['SUCCESS'] = bool(np.all(list(verify.values())))
            if verify['SUCCESS'] is False:
                success = False

            verifyStats[ampName] = verify

        # Loop over amps to make a detector summary.
        verifyDetStats = {'GAIN_FROM_FLAT_PAIR': [], 'ISR_NOISE': []}
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
