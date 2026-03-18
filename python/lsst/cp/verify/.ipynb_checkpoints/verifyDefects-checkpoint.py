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
import scipy.stats
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
from lsst.ip.isr.isrFunctions import countMaskedPixels
import pandas as pd

from .verifyStats import (
    CpVerifyStatsConfig,
    CpVerifyStatsTask,
    CpVerifyStatsConnections,
)


__all__ = ["CpVerifyDefectsConfig", "CpVerifyDefectsTask"]


class CpVerifyDefectsConnections(
    CpVerifyStatsConnections, dimensions={"instrument", "exposure", "detector"}
):
    inputExp = cT.Input(
        name="icExp",
        doc="Input exposure to calculate statistics for.",
        storageClass="Exposure",
        dimensions=["instrument", "exposure", "detector"],
    )
    uncorrectedExp = cT.Input(
        name="uncorrectedExp",
        doc="Uncorrected input exposure to calculate statistics for.",
        storageClass="Exposure",
        dimensions=["instrument", "exposure", "detector"],
    )
    camera = cT.PrerequisiteInput(
        name="camera",
        storageClass="Camera",
        doc="Input camera.",
        dimensions=[
            "instrument",
        ],
        isCalibration=True,
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


class CpVerifyDefectsConfig(
    CpVerifyStatsConfig, pipelineConnections=CpVerifyDefectsConnections
):
    """Inherits from base CpVerifyStatsConfig."""

    def setDefaults(self):
        super().setDefaults()
        self.stageName = "DEFECTS"

        self.maskNameList = ["BAD"]  # noqa F821

        self.imageStatKeywords = {
            "DEFECT_PIXELS": "NMASKED",  # noqa F821
            "OUTLIERS": "NCLIPPED",
            "MEDIAN": "MEDIAN",
            "STDEV": "STDEVCLIP",
            "MIN": "MIN",
            "MAX": "MAX",
        }
        self.unmaskedImageStatKeywords = {
            "UNMASKED_MIN": "MIN",  # noqa F821
            "UNMASKED_MAX": "MAX",
            "UNMASKED_STDEV": "STDEVCLIP",
            "UNMASKED_OUTLIERS": "NCLIPPED",
        }
        self.uncorrectedImageStatKeywords = {
            "UNC_DEFECT_PIXELS": "NMASKED",  # noqa F821
            "UNC_OUTLIERS": "NCLIPPED",
            "UNC_MEDIAN": "MEDIAN",
            "UNC_STDEV": "STDEVCLIP",
            "UNC_MIN": "MIN",
            "UNC_MAX": "MAX",
        }
        # These config options need to have a key/value pair
        # to run verification analysis, but the contents of
        # that pair are not used.
        self.catalogStatKeywords = {"empty": "dictionary"}
        self.detectorStatKeywords = {"empty": "dictionary"}


class CpVerifyDefectsTask(CpVerifyStatsTask):
    """Defects verification sub-class, implementing the verify method.

    This also applies additional image processing statistics.
    """

    ConfigClass = CpVerifyDefectsConfig
    _DefaultName = "cpVerifyDefects"

    def catalogStatistics(self, exposure, catalog, uncorrectedCatalog, statControl):
        """Measure the catalog statistics.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure to measure.
        catalog : `lsst.afw.table.Table`
            The catalog to measure.
        uncorrectedCatalog : `lsst.afw.table.Table`
            The uncorrected catalog to measure.
        statControl : `lsst.afw.math.StatisticsControl`
            Statistics control object with parameters defined by
            the config.

        Returns
        -------
        outputStatistics : `dict` [`str`, `dict` [`str`, scalar]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the statistics measured and their values.

        Notes
        -----
        Number of detections test: running with defects would have fewer
        detections
        """
        outputStatistics = {}

        # Number of detections test
        outputStatistics["NUM_OBJECTS_BEFORE"] = len(uncorrectedCatalog)
        outputStatistics["NUM_OBJECTS_AFTER"] = len(catalog)

        return outputStatistics

    def detectorStatistics(
        self, statisticsDict, statControl, exposure=None, uncorrectedExposure=None
    ):
        """Measure the detector statistics.

        Parameters
        ----------
        statisticsDict : `dict` [`str`, scalar]
            Dictionary with detector tests.
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
            A dictionary containing statistics measured and their values.

        Notes
        -----
        Number of cosmic rays test: If there are defects in our data that
        we didn't properly identify and cover, they might appear similar
        to cosmic rays because they have sharp edges compared to the point
        spread function (PSF). When we process the data, if these defects
        aren't marked in our defect mask, the software might mistakenly think
        they are cosmic rays and try to remove them. However, if we've already
        included these defects in the defect mask, the software won't treat
        them as cosmic rays, so we'll have fewer pixels that are falsely
        identified and removed as cosmic rays when we compare two sets of
        data reductions.
        """
        outputStatistics = {}
        # Cosmic Rays test: Count number of cosmic rays before
        # and after masking with defects
        nCosmicsBefore = countMaskedPixels(uncorrectedExposure, ["CR"])
        nCosmicsAfter = countMaskedPixels(exposure, ["CR"])

        outputStatistics["NUM_COSMICS_BEFORE"] = nCosmicsBefore
        outputStatistics["NUM_COSMICS_AFTER"] = nCosmicsAfter

        return outputStatistics

    def imageStatistics(self, exposure, uncorrectedExposure, statControl):
        """Measure additional defect statistics.

        This calls the parent class method first, then adds additional
        measurements.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure containing the ISR processed data to measure.
        statControl : `lsst.afw.math.StatControl`
            Statistics control object with parameters defined by
            the config.

        Returns
        -------
        outputStatistics : `dict` [`str`, `dict` [`str`, scalar]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the statistics measured and their values.

        Notes
        -----
        Number of cosmic rays test: If there are defects in our data that
        we didn't properly identify and cover, they might appear similar
        to cosmic rays because they have sharp edges compared to the point
        spread function (PSF). When we process the data, if these defects
        aren't marked in our defect mask, the software might mistakenly think
        they are cosmic rays and try to remove them. However, if we've already
        included these defects in the defect mask, the software won't treat
        them as cosmic rays, so we'll have fewer pixels that are falsely
        identified and removed as cosmic rays when we compare two sets of
        data reductions.
        """
        outputStatistics = super().imageStatistics(
            exposure, uncorrectedExposure, statControl
        )

        # Is this a useful test?  It saves having to do chi^2 fits,
        # which are going to be biased by the bulk of points.
        for amp in exposure.getDetector():
            ampName = amp.getName()
            ampExp = exposure.Factory(exposure, amp.getBBox())

            normImage = ampExp.getImage()
            normArray = normImage.getArray()

            normArray -= outputStatistics[ampName]["MEDIAN"]
            normArray /= outputStatistics[ampName]["STDEV"]

            probability = scipy.stats.norm.pdf(normArray)
            outliers = np.where(probability < 1.0 / probability.size, 1.0, 0.0)
            outputStatistics[ampName]["STAT_OUTLIERS"] = int(np.sum(outliers))

            # Get fraction of defects per amp
            ampSize = amp.getBBox().height * amp.getBBox().width
            outputStatistics[ampName]["FRAC"] = (
                outputStatistics[ampName]["DEFECT_PIXELS"] / ampSize
            )

        return outputStatistics

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
            A boolean indicating if all tests have passed.
        """
        # Amplifier statistics
        ampStats = statisticsDict["AMP"]
        verifyStats = {}
        successAmp = True
        for ampName, stats in ampStats.items():
            verify = {}

            # These are not defined in DMTN-101 yet.
            verify["OUTLIERS"] = bool(stats["UNMASKED_OUTLIERS"] >= stats["OUTLIERS"])
            verify["STDEV"] = bool(stats["UNMASKED_STDEV"] >= stats["STDEV"])
            verify["MIN"] = bool(stats["UNMASKED_MIN"] <= stats["MIN"])
            verify["MAX"] = bool(stats["UNMASKED_MAX"] >= stats["MAX"])

            # This test is bad, and should be made not bad.
            verify["PROB_TEST"] = bool(stats["STAT_OUTLIERS"] == stats["DEFECT_PIXELS"])

            verify["SUCCESS"] = bool(np.all(list(verify.values())))
            if verify["SUCCESS"] is False:
                successAmp = False

            verifyStats[ampName] = verify

        success = successAmp
        return {"AMP": verifyStats}, success

    def repackStats(self, statisticsDict, dimensions):
        # docstring inherited
        rows = {}
        rowList = []
        matrixRowList = None

        if self.config.useIsrStatistics:
            mjd = statisticsDict["ISR"]["MJD"]
        else:
            mjd = np.nan

        rowBase = {
            "instrument": dimensions["instrument"],
            "exposure": dimensions[
                "exposure"
            ],  # ensure an exposure dimension for downstream.
            "detector": dimensions["detector"],
            "mjd": mjd,
        }

        # AMP results
        for ampName, stats in statisticsDict["AMP"].items():
            rows[ampName] = {}
            rows[ampName].update(rowBase)
            rows[ampName]["amplifier"] = ampName
            for key, value in stats.items():
                rows[ampName][f"{self.config.stageName}_{key}"] = value

        # ISR results
        nBadColumns = np.nan
        if self.config.useIsrStatistics and "ISR" in statisticsDict:
            for ampName, stats in statisticsDict["ISR"]["CALIBDIST"].items():
                if ampName == "detector":
                    nBadColumns = stats[ampName].get(
                        "LSST CALIB DEFECTS N_BAD_COLUMNS", np.nan
                    )
                elif ampName in stats.keys():
                    key = f"LSST CALIB DEFECTS {ampName} N_HOT"
                    rows[ampName][f"{self.config.stageName}_N_HOT"] = stats[
                        ampName
                    ].get(key, np.nan)
                    key = f"LSST CALIB DEFECTS {ampName} N_COLD"
                    rows[ampName][f"{self.config.stageName}_N_COLD"] = stats[
                        ampName
                    ].get(key, np.nan)

        # DET results
        rows["detector"] = rowBase
        rows["detector"][f"{self.config.stageName}_N_BAD_COLUMNS"] = nBadColumns

        for ampName, stats in rows.items():
            rowList.append(stats)

        return rowList, matrixRowList


# Below is in development
class MeasureDefectsStabilityConnections(
    pipeBase.PipelineTaskConnections, dimensions=("instrument", "detector")
):
    referenceDefects = cT.PrerequisiteInput(
        name="referenceDefects",
        doc="Reference defects.",
        storageClass="Defects",
        dimensions=(
            "instrument",
            "detector",
        ),
        multiple=True,
    )
    productionDefects = cT.PrerequisiteInput(
        name="defectsProduction",
        doc="Measured defects to be compared against reference.",
        storageClass="Defects",
        dimensions=(
            "instrument",
            "detector",
        ),
        multiple=True,
    )
    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera associated with these defects.",
        storageClass="Camera",
        dimensions=("instrument",),
    )
    defectsStability = cT.Output(
        name="defectsStability",
        doc="Dataframe storing the variation of defects.",
        storageClass="DataFrame",
        dimensions=("instrument"),  # TODO: Verify that this dimension is correct
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)


class MeasureDefectsStabilityTaskConfig(
    pipeBase.PipelineTaskConfig, pipelineConnections=MeasureDefectsStabilityConnections
):
    """
    Configuration for measuring stability 
    of defects from individual defect masks.
    """

    connections = pipelineConnections

    refLabel = pexConfig.Field(
        dtype=str,
        doc="Label to apply to the reference defects.",
        optional=False,
    )

    prodLabel = pexConfig.Field(
        dtype=str,
        doc="Label to apply to the production defects.",
        optional=False,
    )

    skipMasks = pexConfig.ListField(
        dtype=str,
        doc="Mask planes to skip.",
        default=(),
        optional=True,
    )

    includeAny = pexConfig.ChoiceField(
        doc="Flag to include the any masks in variability computation.",
        dtype=bool,
        allowed={
            True: "Any mask variability is included in the computation.",
            False: "Any mask variability is excluded in the computation.",
        },
        default=True,
        optional=True,
    )

    def setDefaults(self):

        # Set the acq string based on the two labels
        self.acqString = ""
        for label in [self.refLabel, self.prodLabel]:
            self.acqString += f"{label}_"

        # Set the det_amp objects based on the camera object
        self.detName_Obj = self.connections.camera.getNameMap()
        self.det_amp = []
        for detectorName, detectorObject in self.detName_Obj.items():
            for amp in detectorObject.getAmplifiers():
                self.det_amp.append(f"{detectorName}_{amp.getName()}")

        """
        Here, setup the default dataframe based on the 
        det_amp object and the camera object
        """
        columns = []
        for plane in self.connections.referenceDefects[0].getMaskPlaneDict().keys():
            for acq in [self.refLabel, self.prodLabel]:
                columns.append(f"{acq}_{plane}")
            columns.append(f"{self.acqString}{plane}")

        # If include any flag is set, add the any masks
        if self.includeAny:
            for acq in [self.refLabel, self.prodLabel]:
                columns.append(f"{acq}_Any")
            columns.append(f"{self.acqString}Any")

        self.resultDataFrame = pd.DataFrame(columns=columns, index=self.det_amp)

    def validate(self):
        super().validate()
