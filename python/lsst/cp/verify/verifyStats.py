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
import math

import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
import lsst.meas.algorithms as measAlg

from lsst.cp.pipe.cpCombine import VignetteExposure
from lsst.pipe.tasks.repair import RepairTask
from .utils import mergeStatDict

__all__ = ['CpVerifyStatsConfig', 'CpVerifyStatsTask']


class CpVerifyStatsConnections(pipeBase.PipelineTaskConnections,
                               dimensions={"instrument", "exposure", "detector"},
                               defaultTemplates={}):
    inputExp = cT.Input(
        name="postISRCCD",
        doc="Input exposure to calculate statistics for.",
        storageClass="Exposure",
        dimensions=["instrument", "exposure", "detector"],
    )
    camera = cT.PrerequisiteInput(
        name="camera",
        storageClass="Camera",
        doc="Input camera.",
        dimensions=["instrument", ],
        isCalibration=True,
    )

    outputStats = cT.Output(
        name="detectorStats",
        doc="Output statistics from cp_verify.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "exposure", "detector"],
    )


class CpVerifyStatsConfig(pipeBase.PipelineTaskConfig,
                          pipelineConnections=CpVerifyStatsConnections):
    """Configuration parameters for CpVerifyTask
    """
    maskNameList = pexConfig.ListField(
        dtype=str,
        doc="Mask list to exclude from statistics calculations.",
        default=['DETECTED', 'BAD', 'NO_DATA'],
    )
    doVignette = pexConfig.Field(
        dtype=bool,
        doc="Mask vignetted regions?",
        default=False,
    )
    doNormalize = pexConfig.Field(
        dtype=bool,
        doc="Normalize by exposure time?",
        default=False,
    )

    doCR = pexConfig.Field(
        dtype=bool,
        doc="Run CR rejection?",
        default=False,
    )
    repair = pexConfig.ConfigurableField(
        target=RepairTask,
        doc="Repair task to use.",
    )
    psfFwhm = pexConfig.Field(
        dtype=float,
        default=3.0,
        doc="Repair PSF FWHM (pixels).",
    )
    psfSize = pexConfig.Field(
        dtype=int,
        default=21,
        doc="Repair PSF size (pixels).",
    )
    crGrow = pexConfig.Field(
        dtype=int,
        default=2,
        doc="Grow radius for CR (pixels).",
    )

    useReadNoise = pexConfig.Field(
        dtype=bool,
        doc="Compare sigma against read noise?",
        default=True,
    )

    numSigmaClip = pexConfig.Field(
        dtype=float,
        doc="Rejection threshold (sigma) for statistics clipping.",
        default=3.0,
    )
    clipMaxIter = pexConfig.Field(
        dtype=int,
        doc="Max number of clipping iterations to apply.",
        default=3,
    )

    imageStatKeywords = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Image statistics to run on amplifier segments.",
        default={},
    )
    unmaskedImageStatKeywords = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Image statistics to run on amplifier segments, ignoring masks.",
        default={},
    )
    crImageStatKeywords = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Image statistics to run on CR cleaned amplifier segments.",
        default={},
    )
    normImageStatKeywords = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Image statistics to run on expTime normalized amplifier segments.",
        default={},
    )
    catalogStatKeywords = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Image statistics to run on expTime normalized amplifier segments.",
        default={},
    )
    detectorStatKeywords = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Detector statistics to create.",
        default={},
    )


class CpVerifyStatsTask(pipeBase.PipelineTask, pipeBase.CmdLineTask):
    """
    """
    ConfigClass = CpVerifyStatsConfig
    _DefaultName = 'cpVerifyStats'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.makeSubtask("repair")

    def run(self, inputExp, camera):
        """Calculate quality statistics and verify they meet the requirements
        for a calibration.

        Parameters
        ----------
        inputExp : `lsst.afw.image.Exposure`
            The ISR processed exposure to me measured.
        camera : `lsst.afw.cameraGeom.Camera`
             The camera geometry for this exposure.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:
            - ``outputStats`` : `dict`
                The output measured statistics.

        Notes
        -----
        The outputStats should have a yaml representation of the form

        AMP:
          Amp1:
            STAT: value
            STAT2: value2
          Amp2:
          Amp3:
        DET:
          STAT: value
          STAT2: value
        CATALOG:
          STAT: value
          STAT2: value
        VERIFY:
          DET:
            TEST: boolean
          CATALOG:
            TEST: boolean
          AMP:
            Amp1:
              TEST: boolean
              TEST2: boolean
            Amp2:
            Amp3:
        SUCCESS: boolean

        """
        outputStats = dict()

        if self.config.doVignette:
            VignetteExposure(inputExp, doUpdateMask=True, maskPlane='BAD',
                             doSetValue=False, log=self.log)

        mask = inputExp.getMask()
        maskVal = mask.getPlaneBitMask(self.config.maskNameList)
        statsControl = afwMath.StatisticsControl(self.config.numSigmaClip,
                                                 self.config.clipMaxIter,
                                                 maskVal)

        # This is wrapped below to check for config lengths, as we can
        # make a number of different image stats.
        outputStats['AMP'] = self.imageStatistics(inputExp, statsControl)
        if len(self.config.catalogStatKeywords):
            outputStats['CATALOG'] = self.catalogStatistics(inputExp, statsControl)
        else:
            outputStats['CATALOG'] = dict()
        if len(self.config.detectorStatKeywords):
            outputStats['DET'] = self.detectorStatistics(outputStats, statsControl)
        else:
            outputStats['DET'] = dict()

        outputStats['VERIFY'], outputStats['SUCCESS'] = self.verify(inputExp, outputStats)

        return pipeBase.Struct(
            outputStats=outputStats,
        )

    @staticmethod
    def _emptyAmpDict(exposure):
        """Construct empty dictionary indexed by amplifier names.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to extract detector from.

        Returns
        -------
        outputStats : `dict` [`str`, `dict`]
            A skeleton statistics dictionary.

        Raises
        ------
        RuntimeError :
            Raised if no detector can be found.
        """
        outputStats = dict()
        detector = exposure.getDetector()
        if detector is None:
            raise RuntimeError("No detector found in exposure!")

        for amp in detector.getAmplifiers():
            outputStats[amp.getName()] = dict()

        return outputStats

    # Image measurement methods.
    def imageStatistics(self, exposure, statControl):
        """Measure all types of image statistics.

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

        outputStats = self._emptyAmpDict(exposure)

        if len(self.config.imageStatKeywords):
            outputStats = mergeStatDict(outputStats, self.amplifierStats(exposure,
                                                                         self.config.imageStatKeywords,
                                                                         statControl))

        if len(self.config.unmaskedImageStatKeywords):
            outputStats = mergeStatDict(outputStats, self.unmaskedImageStats(exposure, statControl))

        if len(self.config.normImageStatKeywords):
            outputStats = mergeStatDict(outputStats, self.normalizedImageStats(exposure, statControl))

        if len(self.config.crImageStatKeywords):
            outputStats = mergeStatDict(outputStats, self.crImageStats(exposure, statControl))

        return outputStats

    def amplifierStats(self, exposure, keywordDict, statsControl):
        """Measure amplifier level statistics on the data read from disk.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure to measure.
        keywordDict : `dict` [`str`, `str`]
            A dictionary of keys to use in the output results, with
            values the string name associated with the
            `lsst.afw.math.statistics.Property` to measure.
        statControl : `lsst.afw.math.StatControl`
            Statistics control object with the parameters defined by
            the config.

        Returns
        -------
        ampStats : `dict` [`str`, `dict` [`str`, scalar]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the statistics measured and their values.
        """

        ampStats = dict()

        # Convert the string names in the keywordDict to the afwMath values.
        # The statisticToRun is then the bitwise-or of that set.
        statisticToRun = 0
        statAccessor = {}
        for k, v in keywordDict.items():
            statValue = afwMath.stringToStatisticsProperty(v)
            statisticToRun |= statValue
            statAccessor[k] = statValue

        # Measure stats on all amplifiers.
        for ampIdx, amp in enumerate(exposure.getDetector()):
            ampName = amp.getName()
            theseStats = dict()
            ampExp = exposure.Factory(exposure, amp.getBBox())
            stats = afwMath.makeStatistics(ampExp.getMaskedImage(), statisticToRun, statsControl)

            for k, v in statAccessor.items():
                theseStats[k] = stats.getValue(v)
            ampStats[ampName] = theseStats

        return ampStats

    def unmaskedImageStats(self, exposure, statControl):
        """Measure amplifier level statistics on the data read from disk,
        ignoring all mask planes.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure to measure.
        statControl : `lsst.afw.math.StatControl`
            Statistics control object with the parameters defined by
            the config.

        Returns
        -------
        outputStats : `dict` [`str`, `dict` [`str`, scalar]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the statistics measured and their values.

        """
        noMaskStatsControl = afwMath.StatisticsControl(self.config.numSigmaClip,
                                                       self.config.clipMaxIter,
                                                       0x0)
        return self.amplifierStats(exposure, self.config.unmaskedImageStatKeywords, noMaskStatsControl)

    def normalizedImageStats(self, exposure, statControl):
        """Measure amplifier level statistics on the data read from disk,
        after dividing by the exposure time.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure to measure.
        statControl : `lsst.afw.math.StatControl`
            Statistics control object with the parameters defined by
            the config.

        Returns
        -------
        outputStats : `dict` [`str`, `dict` [`str`, scalar]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the statistics measured and their values.

        """
        scaledExposure = exposure.clone()
        exposureTime = scaledExposure.getInfo().getVisitInfo().getExposureTime()
        mi = scaledExposure.getMaskedImage()
        mi /= exposureTime

        return self.amplifierStats(scaledExposure, self.config.normImageStatKeywords, statControl)

    def crImageStats(self, exposure, statControl):
        """Measure amplifier level statistics on the data read from disk,
        after running cosmic ray rejection.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure to measure.
        statControl : `lsst.afw.math.StatControl`
            Statistics control object with the parameters defined by
            the config.

        Returns
        -------
        outputStats : `dict` [`str`, `dict` [`str`, scalar]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the statistics measured and their values.

        """
        crRejectedExp = exposure.clone()
        psf = measAlg.SingleGaussianPsf(self.config.psfSize,
                                        self.config.psfSize,
                                        self.config.psfFwhm/(2*math.sqrt(2*math.log(2))))
        crRejectedExp.setPsf(psf)
        self.repair.run(crRejectedExp, keepCRs=False)
        if self.config.crGrow > 0:
            crMask = crRejectedExp.getMaskedImage().getMask().getPlaneBitMask("CR")
            spans = afwGeom.SpanSet.fromMask(crRejectedExp.mask, crMask)
            spans = spans.dilated(self.config.crGrow)
            spans.setMask(crRejectedExp.mask, crMask)

        return self.amplifierStats(crRejectedExp, self.config.crImageStatKeywords, statControl)

    # Methods that need to be implemented by the calibration-level subclasses.
    def catalogStatistics(self, exposure, statControl):
        """Calculate statistics from a catalog.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure to measure.
        statControl : `lsst.afw.math.StatControl`
            Statistics control object with the parameters defined by
            the config.

        Returns
        -------
        outputStats : `dict` [`str`, `dict` [`str`, scalar]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the statistics measured and their values.
        """
        raise NotImplementedError("Subclasses must implement verification criteria.")

    def detectorStatistics(self, statisticsDictionary, statControl):
        """Calculate detector level statistics based on the existing measurements.

        Parameters
        ----------
        statisticsDictionary : `dict` [`str`, `dict` [`str`, scalar]],
            Dictionary of measured statistics.

        Returns
        -------
        outputStats : `dict` [`str`, scalar]
            A dictionary of the statistics measured and their values.

        Raises
        ------
        NotImplementedError :
            This method must be implemented by the calibration-type
            subclass.
        """
        raise NotImplementedError("Subclasses must implement verification criteria.")

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

        Raises
        ------
        NotImplementedError :
            This method must be implemented by the calibration-type
            subclass.
        """
        raise NotImplementedError("Subclasses must implement verification criteria.")
