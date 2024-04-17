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

from astropy.table import Table

import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pex.exceptions as pexException
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
import lsst.meas.algorithms as measAlg

from lsst.ip.isr.vignette import maskVignettedRegion
from lsst.pipe.tasks.repair import RepairTask
from .utils import mergeStatDict


__all__ = ["CpVerifyStatsConfig", "CpVerifyStatsTask"]


class CpVerifyStatsConnections(
    pipeBase.PipelineTaskConnections,
    dimensions={"instrument", "exposure", "detector"},
    defaultTemplates={},
):
    inputExp = cT.Input(
        name="postISRCCD",
        doc="Input exposure to calculate statistics for.",
        storageClass="Exposure",
        dimensions=["instrument", "exposure", "detector"],
    )
    uncorrectedExp = cT.Input(
        name="uncorrectedExp",
        doc="Uncorrected input exposure to calculate statistics for.",
        storageClass="ExposureF",
        dimensions=["instrument", "visit", "detector"],
    )
    taskMetadata = cT.Input(
        name="isrTask_metadata",
        doc="Input task metadata to extract statistics from.",
        storageClass="TaskMetadata",
        dimensions=["instrument", "exposure", "detector"],
    )
    inputCatalog = cT.Input(
        name="src",
        doc="Input catalog to calculate statistics for.",
        storageClass="SourceCatalog",
        dimensions=["instrument", "visit", "detector"],
    )
    uncorrectedCatalog = cT.Input(
        name="uncorrectedSrc",
        doc="Input catalog without correction applied.",
        storageClass="SourceCatalog",
        dimensions=["instrument", "visit", "detector"],
    )
    camera = cT.PrerequisiteInput(
        name="camera",
        storageClass="Camera",
        doc="Input camera.",
        dimensions=["instrument", ],
        isCalibration=True,
    )
    isrStatistics = cT.Input(
        name="isrStatistics",
        storageClass="StructuredDataDict",
        doc="Pre-calculated statistics from IsrTask.",
        dimensions=["instrument", "exposure", "detector"],
    )

    outputStats = cT.Output(
        name="detectorStats",
        doc="Output statistics from cp_verify.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "exposure", "detector"],
    )
    outputResults = cT.Output(
        name="detectorResults",
        doc="Output results from cp_verify.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "exposure", "detector"],
    )
    outputMatrix = cT.Output(
        name="detectorMatrix",
        doc="Output matrix results from cp_verify.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "exposure", "detector"],
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if len(config.metadataStatKeywords) < 1:
            self.inputs.discard("taskMetadata")

        if len(config.catalogStatKeywords) < 1:
            self.inputs.discard("inputCatalog")
            self.inputs.discard("uncorrectedCatalog")

        if len(config.uncorrectedImageStatKeywords) < 1:
            self.inputs.discard("uncorrectedExp")

        if config.useIsrStatistics is not True:
            self.inputs.discard("isrStatistics")

        if not config.hasMatrixCatalog:
            self.outputs.discard("outputMatrix")


class CpVerifyStatsConfig(
    pipeBase.PipelineTaskConfig, pipelineConnections=CpVerifyStatsConnections
):
    """Configuration parameters for CpVerifyStatsTask."""

    maskNameList = pexConfig.ListField(
        dtype=str,
        doc="Mask list to exclude from statistics calculations.",
        default=["DETECTED", "BAD", "NO_DATA"],
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

    # Cosmic ray handling options.
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
        doc="Repair PSF bounding-box size (pixels).",
    )
    crGrow = pexConfig.Field(
        dtype=int,
        default=0,
        doc="Grow radius for CR (pixels).",
    )

    # Statistics options.
    useReadNoise = pexConfig.Field(
        dtype=bool,
        doc="Compare sigma against read noise?",
        default=True,
    )
    numSigmaClip = pexConfig.Field(
        dtype=float,
        doc="Rejection threshold (sigma) for statistics clipping.",
        default=5.0,
    )
    clipMaxIter = pexConfig.Field(
        dtype=int,
        doc="Max number of clipping iterations to apply.",
        default=3,
    )

    # Keywords and statistics to measure from different sources.
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
    uncorrectedImageStatKeywords = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Uncorrected image statistics to run on amplifier segments.",
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
    metadataStatKeywords = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Statistics to measure from the metadata of the exposure.",
        default={},
    )
    catalogStatKeywords = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Statistics to measure from source catalogs of objects in the exposure.",
        default={},
    )
    detectorStatKeywords = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Statistics to create for the full detector from the per-amplifier measurements.",
        default={},
    )

    stageName = pexConfig.Field(
        dtype=str,
        doc="Stage name to use for table columns.",
        default="NOSTAGE",
    )
    useIsrStatistics = pexConfig.Field(
        dtype=bool,
        doc="Use statistics calculated by IsrTask?",
        default=False,
    )
    hasMatrixCatalog = pexConfig.Field(
        dtype=bool,
        doc="Will a matrix table of results be made?",
        default=False,
    )
    expectedDistributionLevels = pexConfig.ListField(
        dtype=float,
        doc="Percentile levels expected in the calibration header.",
        default=[0, 5, 16, 50, 84, 95, 100],
    )


class CpVerifyStatsTask(pipeBase.PipelineTask):
    """Main statistic measurement and validation class.

    This operates on a single (exposure, detector) pair, and is
    designed to be subclassed so specific calibrations can apply their
    own validation methods.
    """

    ConfigClass = CpVerifyStatsConfig
    _DefaultName = "cpVerifyStats"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.makeSubtask("repair")

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        # Pass the full dataId, as we want to retain filter info.
        inputs["dimensions"] = dict(inputRefs.inputExp.dataId.mapping)
        print("CZW", inputs["dimensions"], dict(inputRefs.inputExp.dataId.mapping))
        # import pdb; pdb.set_trace()
        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def run(
        self,
        inputExp,
        camera,
        isrStatistics=None,
        uncorrectedExp=None,
        taskMetadata=None,
        inputCatalog=None,
        uncorrectedCatalog=None,
        dimensions=None,
    ):
        """Calculate quality statistics and verify they meet the requirements
        for a calibration.

        Parameters
        ----------
        inputExp : `lsst.afw.image.Exposure`
            The ISR processed exposure to be measured.
        camera : `lsst.afw.cameraGeom.Camera`
            The camera geometry for ``inputExp``.
        uncorrectedExp : `lsst.afw.image.Exposure`
            The alternate exposure to measure.
        taskMetadata : `lsst.pipe.base.TaskMetadata`, optional
            Task metadata containing additional statistics.
        inputCatalog : `lsst.afw.image.Table`
            The source catalog to measure.
        uncorrectedCatalog : `lsst.afw.image.Table`
            The alternate source catalog to measure.
        dimensions : `dict`
            Dictionary of input dictionary.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:
            - ``outputStats`` : `dict`
                The output measured statistics.
        """
        outputStats = {}

        if self.config.doVignette:
            polygon = inputExp.getInfo().getValidPolygon()
            maskVignettedRegion(
                inputExp, polygon, maskPlane="NO_DATA", vignetteValue=None, log=self.log
            )

        mask = inputExp.getMask()
        maskVal = mask.getPlaneBitMask(self.config.maskNameList)
        statControl = afwMath.StatisticsControl(
            self.config.numSigmaClip, self.config.clipMaxIter, maskVal
        )

        # This is wrapped below to check for config lengths, as we can
        # make a number of different image stats.
        outputStats["AMP"] = self.imageStatistics(inputExp, uncorrectedExp, statControl)

        if len(self.config.metadataStatKeywords):
            # These are also defined on a amp-by-amp basis.
            outputStats["METADATA"] = self.metadataStatistics(inputExp, taskMetadata)
        else:
            outputStats["METADATA"] = {}

        if len(self.config.catalogStatKeywords):
            outputStats["CATALOG"] = self.catalogStatistics(
                inputExp, inputCatalog, uncorrectedCatalog, statControl
            )
        else:
            outputStats["CATALOG"] = {}
        if len(self.config.detectorStatKeywords):
            outputStats["DET"] = self.detectorStatistics(
                outputStats, statControl, inputExp, uncorrectedExp
            )
        else:
            outputStats["DET"] = {}

        if self.config.useIsrStatistics:
            outputStats["ISR"] = isrStatistics

        outputStats["VERIFY"], outputStats["SUCCESS"] = self.verify(
            inputExp, outputStats
        )

        outputResults, outputMatrix = self.repackStats(outputStats, dimensions)

        return pipeBase.Struct(
            outputStats=outputStats,
            outputResults=Table(outputResults),
            outputMatrix=Table(outputMatrix),
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
        outputStatistics : `dict` [`str`, `dict`]
            A skeleton statistics dictionary.

        Raises
        ------
        RuntimeError :
            Raised if no detector can be found.
        """
        outputStatistics = {}
        detector = exposure.getDetector()
        if detector is None:
            raise RuntimeError("No detector found in exposure!")

        for amp in detector.getAmplifiers():
            outputStatistics[amp.getName()] = {}

        return outputStatistics

    # Image measurement methods.
    def imageStatistics(self, exposure, uncorrectedExposure, statControl):
        """Measure image statistics for a number of simple image
        modifications.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure containing the ISR processed data to measure.
        uncorrectedExposure: `lsst.afw.image.Exposure`
            Uncorrected exposure containing the ISR processed data to measure.
        statControl : `lsst.afw.math.StatisticsControl`
            Statistics control object with parameters defined by
            the config.

        Returns
        -------
        outputStatistics : `dict` [`str`, `dict` [`str`, scalar]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the statistics measured and their values.

        """
        outputStatistics = self._emptyAmpDict(exposure)

        if len(self.config.imageStatKeywords):
            outputStatistics = mergeStatDict(
                outputStatistics,
                self.amplifierStats(
                    exposure, self.config.imageStatKeywords, statControl
                ),
            )
        if len(self.config.uncorrectedImageStatKeywords):
            outputStatistics = mergeStatDict(
                outputStatistics,
                self.amplifierStats(
                    uncorrectedExposure,
                    self.config.uncorrectedImageStatKeywords,
                    statControl,
                ),
            )
        if len(self.config.unmaskedImageStatKeywords):
            outputStatistics = mergeStatDict(
                outputStatistics, self.unmaskedImageStats(exposure)
            )

        if len(self.config.normImageStatKeywords):
            outputStatistics = mergeStatDict(
                outputStatistics, self.normalizedImageStats(exposure, statControl)
            )

        if len(self.config.crImageStatKeywords):
            outputStatistics = mergeStatDict(
                outputStatistics, self.crImageStats(exposure, statControl)
            )

        return outputStatistics

    @staticmethod
    def _configHelper(keywordDict):
        """Helper to convert keyword dictionary to stat value.

        Convert the string names in the keywordDict to the afwMath values.
        The statisticToRun is then the bitwise-or of that set.

        Parameters
        ----------
        keywordDict : `dict` [`str`, `str`]
            A dictionary of keys to use in the output results, with
            values the string name associated with the
            `lsst.afw.math.statistics.Property` to measure.

        Returns
        -------
        statisticToRun : `int`
            The merged `lsst.afw.math` statistics property.
        statAccessor : `dict` [`str`, `int`]
            Dictionary containing statistics property indexed by name.
        """
        statisticToRun = 0
        statAccessor = {}
        for k, v in keywordDict.items():
            statValue = afwMath.stringToStatisticsProperty(v)
            statisticToRun |= statValue
            statAccessor[k] = statValue

        return statisticToRun, statAccessor

    def metadataStatistics(self, exposure, taskMetadata):
        """Extract task metadata information for verification.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure to measure.
        taskMetadata : `lsst.pipe.base.TaskMetadata`
            The metadata to extract values from.

        Returns
        -------
        ampStats : `dict` [`str`, `dict` [`str`, scalar]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the statistics measured and their values.
        """
        metadataStats = {}
        keywordDict = self.config.metadataStatKeywords

        if taskMetadata:
            for key, value in keywordDict.items():
                if value == "AMP":
                    metadataStats[key] = {}
                    for ampIdx, amp in enumerate(exposure.getDetector()):
                        ampName = amp.getName()
                        expectedKey = f"{key} {ampName}"
                        metadataStats[key][ampName] = None
                        for name in taskMetadata:
                            if expectedKey in taskMetadata[name]:
                                metadataStats[key][ampName] = taskMetadata[name][
                                    expectedKey
                                ]
                else:
                    # Assume it's detector-wide.
                    expectedKey = key
                    for name in taskMetadata:
                        if expectedKey in taskMetadata[name]:
                            metadataStats[key] = taskMetadata[name][expectedKey]
        return metadataStats

    def amplifierStats(self, exposure, keywordDict, statControl, failAll=False):
        """Measure amplifier level statistics from the exposure.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure to measure.
        keywordDict : `dict` [`str`, `str`]
            A dictionary of keys to use in the output results, with
            values the string name associated with the
            `lsst.afw.math.statistics.Property` to measure.
        statControl : `lsst.afw.math.StatisticsControl`
            Statistics control object with parameters defined by
            the config.
        failAll : `bool`, optional
            If True, all tests will be set as failed.

        Returns
        -------
        ampStats : `dict` [`str`, `dict` [`str`, scalar]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the statistics measured and their values.
        """
        ampStats = {}
        statisticToRun, statAccessor = self._configHelper(keywordDict)
        # Measure stats on all amplifiers.
        for ampIdx, amp in enumerate(exposure.getDetector()):
            ampName = amp.getName()
            theseStats = {}
            ampExp = exposure.Factory(exposure, amp.getBBox())
            stats = afwMath.makeStatistics(
                ampExp.getMaskedImage(), statisticToRun, statControl
            )

            for k, v in statAccessor.items():
                theseStats[k] = stats.getValue(v)

            if failAll:
                theseStats["FORCE_FAILURE"] = failAll
            ampStats[ampName] = theseStats

        return ampStats

    def unmaskedImageStats(self, exposure):
        """Measure amplifier level statistics on the exposure, including all
        pixels in the exposure, regardless of any mask planes set.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure to measure.

        Returns
        -------
        outputStatistics : `dict` [`str`, `dict` [`str`, scalar]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the statistics measured and their values.
        """
        noMaskStatsControl = afwMath.StatisticsControl(
            self.config.numSigmaClip, self.config.clipMaxIter, 0x0
        )
        return self.amplifierStats(
            exposure, self.config.unmaskedImageStatKeywords, noMaskStatsControl
        )

    def normalizedImageStats(self, exposure, statControl):
        """Measure amplifier level statistics on the exposure after dividing
        by the exposure time.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure to measure.
        statControl : `lsst.afw.math.StatisticsControl`
            Statistics control object with parameters defined by
            the config.

        Returns
        -------
        outputStatistics : `dict` [`str`, `dict` [`str`, scalar]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the statistics measured and their values.

        Raises
        ------
        RuntimeError :
            Raised if the exposure time cannot be used for normalization.
        """
        scaledExposure = exposure.clone()
        exposureTime = scaledExposure.getInfo().getVisitInfo().getExposureTime()
        if exposureTime <= 0:
            raise RuntimeError(f"Invalid exposureTime {exposureTime}.")
        mi = scaledExposure.getMaskedImage()
        mi /= exposureTime

        return self.amplifierStats(
            scaledExposure, self.config.normImageStatKeywords, statControl
        )

    def crImageStats(self, exposure, statControl):
        """Measure amplifier level statistics on the exposure,
        after running cosmic ray rejection.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure to measure.
        statControl : `lsst.afw.math.StatisticsControl`
            Statistics control object with parameters defined by
            the config.

        Returns
        -------
        outputStatistics : `dict` [`str`, `dict` [`str`, scalar]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the statistics measured and their values.

        """
        crRejectedExp = exposure.clone()
        psf = measAlg.SingleGaussianPsf(
            self.config.psfSize,
            self.config.psfSize,
            self.config.psfFwhm / (2 * np.sqrt(2 * np.log(2))),
        )
        crRejectedExp.setPsf(psf)
        try:
            self.repair.run(crRejectedExp, keepCRs=False)
            failAll = False
        except pexException.LengthError:
            self.log.warning(
                "Failure masking cosmic rays (too many found).  Continuing."
            )
            failAll = True

        if self.config.crGrow > 0:
            crMask = crRejectedExp.getMaskedImage().getMask().getPlaneBitMask("CR")
            spans = afwGeom.SpanSet.fromMask(crRejectedExp.mask, crMask)
            spans = spans.dilated(self.config.crGrow)
            spans = spans.clippedTo(crRejectedExp.getBBox())
            spans.setMask(crRejectedExp.mask, crMask)

        return self.amplifierStats(
            crRejectedExp, self.config.crImageStatKeywords, statControl, failAll=failAll
        )

    # Methods that need to be implemented by the calibration-level subclasses.
    def catalogStatistics(self, exposure, catalog, uncorrectedCatalog, statControl):
        """Calculate statistics from a catalog.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure to measure.
        catalog : `lsst.afw.table.Table`
            The catalog to measure.
        uncorrectedCatalog : `lsst.afw.table.Table`
            The alternate catalog to measure.
        statControl : `lsst.afw.math.StatisticsControl`
            Statistics control object with parameters defined by
            the config.

        Returns
        -------
        outputStatistics : `dict` [`str`, `dict` [`str`, scalar]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the statistics measured and their values.
        """
        raise NotImplementedError(
            "Subclasses must implement catalog statistics method."
        )

    def detectorStatistics(
        self, statisticsDict, statControl, exposure=None, uncorrectedExposure=None
    ):
        """Calculate detector level statistics based on the existing
        per-amplifier measurements.

        Parameters
        ----------
        statisticsDict : `dict` [`str`, scalar]
            Dictionary of measured statistics.  The inner dictionary
            should have keys that are statistic names (`str`) with
            values that are some sort of scalar (`int` or `float` are
            the mostly likely types).
        statControl : `lsst.afw.math.StatControl`
            Statistics control object with parameters defined by
            the config.
        exposure : `lsst.afw.image.Exposure`, optional
            Exposure containing the ISR-processed data to measure.
        uncorrectedExposure : `lsst.afw.image.Exposure`, optional
            uncorrected esposure (no defects) containing the
            ISR-processed data to measure.

        Returns
        -------
        outputStatistics : `dict` [`str`, scalar]
            A dictionary of the statistics measured and their values.

        Raises
        ------
        NotImplementedError :
            This method must be implemented by the calibration-type
            subclass.
        """
        raise NotImplementedError(
            "Subclasses must implement detector statistics method."
        )

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
            A boolean indicating whether all tests have passed.

        Raises
        ------
        NotImplementedError :
            This method must be implemented by the calibration-type
            subclass.
        """
        raise NotImplementedError("Subclasses must implement verification criteria.")

    def repackStats(self, statisticsDict, dimensions):
        """Repack information into flat tables.

        This method may be redefined in subclasses.  This default
        version will repack simple amp-level statistics and
        verification results.

        Parameters
        ----------
        statisticsDictionary : `dict` [`str`, `dict` [`str`, scalar]],
            Dictionary of measured statistics.  The inner dictionary
            should have keys that are statistic names (`str`) with
            values that are some sort of scalar (`int` or `float` are
            the mostly likely types).
        dimensions : `dict`
            The dictionary of dimensions values for this data, to be
            included in the output results.

        Returns
        -------
        outputResults : `list` [`dict`]
            A list of rows to add to the output table.
        outputMatrix : `list` [`dict`]
            A list of rows to add to the output matrix.
        """
        rows = {}
        rowList = []
        matrixRowList = None

        if self.config.useIsrStatistics:
            mjd = statisticsDict["ISR"]["MJD"]
        else:
            mjd = np.nan

        rowBase = {
            "instrument": dimensions["instrument"],
            "detector": dimensions["detector"],
            "mjd": mjd,
        }

        # AMP results:
        for ampName, stats in statisticsDict["AMP"].items():
            rows[ampName] = {}
            rows[ampName].update(rowBase)
            rows[ampName]["amplifier"] = ampName
            for key, value in stats.items():
                rows[ampName][f"{self.config.stageName}_{key}"] = value

        # VERIFY results
        if "AMP" in statisticsDict["VERIFY"]:
            for ampName, stats in statisticsDict["VERIFY"]["AMP"].items():
                for key, value in stats.items():
                    rows[ampName][f"{self.config.stageName}_VERIFY_{key}"] = value

        # pack final list
        for ampName, stats in rows.items():
            rowList.append(stats)

        return rowList, matrixRowList
