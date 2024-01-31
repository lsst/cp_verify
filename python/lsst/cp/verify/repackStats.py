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

import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
import lsst.pex.config as pexConfig

__all__ = [
    "CpVerifyRepackInstrumentConnections",
    "CpVerifyRepackPhysicalFilterConnections",
    "CpVerifyRepackInstrumentConfig",
    "CpVerifyRepackPhysicalFilterConfig",
    "CpVerifyRepackBiasTask",
    "CpVerifyRepackDarkTask",
    "CpVerifyRepackFlatTask",
    "CpVerifyRepackDefectTask",
    "CpVerifyRepackPtcTask",
    "CpVerifyRepackBfkTask",
    "CpVerifyRepackCtiTask",
]


class CpVerifyRepackInstrumentConnections(pipeBase.PipelineTaskConnections,
                                          dimensions={"instrument"},
                                          defaultTemplate={}):
    """Connections class for calibration statistics with only instrument
    dimension.
    """
    detectorStats = cT.Input(
        name="detectorStats",
        doc="Input detector statistics.",
        storageClass="StructuredDataDict",
        dimensions={"instrument", "exposure", "detector"},
        multiple=True,
    )
    exposureStats = cT.Input(
        name="exposureStats",
        doc="Input exposure statistics.",
        storageClass="StructuredDataDict",
        dimensions={"instrument", "exposure"},
        multiple=True,
    )
    runStats = cT.Input(
        name="runStats",
        doc="Input Run statistics.",
        storageClass="StructuredDataDict",
        dimensions={"instrument", },
        multiple=True,
    )

    outputCatalog = cT.Output(
        name="cpvCatalog",
        doc="Output merged catalog.",
        storageClass="ArrowAstropy",
        dimensions={"instrument", },
    )
    matrixCatalog = cT.Output(
        name="cpvMatrix",
        doc="Output matrix catalog.",
        storageClass="ArrowAstropy",
        dimensions={"instrument", },
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not config.hasMatrixCatalog:
            self.outputs.remove("matrixCatalog")


class CpVerifyRepackInstrumentConfig(pipeBase.PipelineTaskConfig,
                                     pipelineConnections=CpVerifyRepackInstrumentConnections):

    expectedDistributionLevels = pexConfig.ListField(
        dtype=float,
        doc="Percentile levels expected in the calibration header.",
        default=[0, 5, 16, 50, 84, 95, 100],
    )
    hasMatrixCatalog = pexConfig.Field(
        dtype=bool,
        doc="Will a matrix catalog be created?",
        default=False,
    )


class CpVerifyRepackTask(pipeBase.PipelineTask):
    """Repack cpVerify statistics for analysis_tools.

    This version is the base for calibrations with summary
    dimensions of instrument only.
    """
    ConfigClass = CpVerifyRepackInstrumentConfig
    _DefaultName = "cpVerifyRepack"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        inputs["detectorDims"] = [dict(exp.dataId.required) for exp in inputRefs.detectorStats]
        inputs["exposureDims"] = [dict(exp.dataId.required) for exp in inputRefs.exposureStats]

        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def run(self, detectorStats, detectorDims, exposureStats, exposureDims, runStats):
        """
        """
        results = self.repack(detectorStats, detectorDims, exposureStats, exposureDims, runStats)
        catalog = Table(results.rowList)

        matrixCatalog = None
        if self.config.hasMatrixCatalog:
            matrixCatalog = Table(results.matrixList)
        return pipeBase.Struct(
            outputCatalog=catalog,
            matrixCatalog=matrixCatalog,
        )

    def repackDetStats(self, detectorStats, detectorDims):
        raise NotImplementedError("Repack needs to be defined by subclasses.")

    def repackExpStats(self, exposureStats, exposureDims):
        # for expStats, expDims in zip(exposureStats, exposureDims):
        raise NotImplementedError("Repack needs to be defined by subclasses.")

    def repackRunStats(self, runStats):
        # for runStats in runStats:
        raise NotImplementedError("Repack needs to be defined by subclasses.")

    def repack(self, detectorStats, detectorDims, exposureStats, exposureDims, runStats):
        return self.repackDetStats(detectorStats, detectorDims)


class CpVerifyRepackBiasTask(CpVerifyRepackTask):
    stageName = "bias"

    def repackDetStats(self, detectorStats, detectorDims):
        rowList = []

        for detStats, detDims in zip(detectorStats, detectorDims):
            row = {}
            instrument = detDims["instrument"]
            exposure = detDims["exposure"]
            detector = detDims["detector"]
            mjd = detStats["ISR"]["MJD"]

            # Get amp stats
            # AMP {ampName} [CR_NOISE MEAN NOISE] value
            for ampName, stats in detStats["AMP"].items():
                row[ampName] = {
                    "instrument": instrument,
                    "exposure": exposure,
                    "mjd": mjd,
                    "detector": detector,
                    "amplifier": ampName,
                    "biasMean": stats["MEAN"],
                    "biasNoise": stats["NOISE"],
                    "biasCrNoise": stats["CR_NOISE"]
                }
            # Get catalog stats  CATALOG
            # Get detector stats DET
            # Get metadata stats
            # METADATA (RESIDUAL STDEV) {ampName} value
            for ampName, value in detStats["METADATA"]["RESIDUAL STDEV"].items():
                row[ampName]["biasReadNoise"] = value

            # Get verify stats
            for ampName, stats in detStats["VERIFY"]["AMP"].items():
                row[ampName]["biasVerifyMean"] = stats["MEAN"]
                row[ampName]["biasVerifyNoise"] = stats["NOISE"]
                row[ampName]["biasVerifyCrNoise"] = stats["CR_NOISE"]
                row[ampName]["biasVerifyReadNoiseConsistent"] = stats["READ_NOISE_CONSISTENT"]

            # Get isr stats
            for ampName, stats in detStats["ISR"]["CALIBDIST"].items():
                for level in self.config.expectedDistributionLevels:
                    key = f"LSST CALIB {self.stageName.upper()} {ampName} DISTRIBUTION {level}-PCT"
                    row[ampName][f"biasDistribution_{level}"] = stats[key]

            projStats = detStats["ISR"]["PROJECTION"]
            for ampName in projStats["AMP_HPROJECTION"].keys():
                row[ampName]["biasSerialProfile"] = np.array(projStats["AMP_HPROJECTION"][ampName])
            for ampName in projStats["AMP_VPROJECTION"].keys():
                row[ampName]["biasParallelProfile"] = np.array(projStats["AMP_VPROJECTION"][ampName])

            shiftStats = detStats["ISR"]["BIASSHIFT"]
            for ampName, stats in shiftStats.items():
                row[ampName]["biasShiftCount"] = len(stats["BIAS_SHIFTS"])
                row[ampName]["biasShiftNoise"] = stats["LOCAL_NOISE"]
            corrStats = detStats["ISR"]["AMPCORR"]

            # Create output table:
            for ampName, stats in row.items():
                rowList.append(stats)

        # We need all rows of biasParallelProfile and biasParallelProfile
        # to be the same length for serialization. Therefore, we pad
        # to the longest length.

        maxSerialLen = 0
        maxParallelLen = 0

        for row in rowList:
            if len(row["biasSerialProfile"]) > maxSerialLen:
                maxSerialLen = len(row["biasSerialProfile"])
            if len(row["biasParallelProfile"]) > maxParallelLen:
                maxParallelLen = len(row["biasParallelProfile"])

        for row in rowList:
            if len(row["biasSerialProfile"]) < maxSerialLen:
                row["biasSerialProfile"] = np.pad(
                    row["biasSerialProfile"],
                    (0, maxSerialLen - len(row["biasSerialProfile"])),
                    constant_values=np.nan,
                )
            if len(row["biasParallelProfile"]) < maxParallelLen:
                row["biasParallelProfile"] = np.pad(
                    row["biasParallelProfile"],
                    (0, maxParallelLen - len(row["biasParallelProfile"])),
                    constant_values=np.nan,
                )

        return pipeBase.Struct(
            rowList=rowList,
            matrixList=corrStats,
        )


class CpVerifyRepackDarkTask(CpVerifyRepackTask):
    stageName = "dark"

    def repackDetStats(self, detectorStats, detectorDims):
        rowList = []

        for detStats, detDims in zip(detectorStats, detectorDims):
            row = {}
            instrument = detDims["instrument"]
            exposure = detDims["exposure"]
            detector = detDims["detector"]

            # Get amp stats
            # AMP {ampName} [CR_NOISE MEAN NOISE] value
            for ampName, stats in detStats["AMP"].items():
                row[ampName] = {
                    "instrument": instrument,
                    "exposure": exposure,
                    "detector": detector,
                    "amplifier": ampName,
                    "darkMean": stats["MEAN"],
                    "darkNoise": stats["NOISE"],
                    "darkCrNoise": stats["CR_NOISE"]
                }
            # Get catalog stats
            # Get detector stats
            # Get metadata stats
            for ampName, value in detStats["METADATA"]["RESIDUAL STDEV"].items():
                row[ampName]["darkReadNoise"] = value

            # Get verify stats
            for ampName, stats in detStats["VERIFY"]["AMP"].items():
                row[ampName]["darkVerifyMean"] = stats["MEAN"]
                row[ampName]["darkVerifyNoise"] = stats["NOISE"]
                row[ampName]["darkVerifyCrNoise"] = stats["CR_NOISE"]
                row[ampName]["darkVerifyReadNoiseConsistent"] = stats["READ_NOISE_CONSISTENT"]
            # Get isr stats
            for ampName, stats in detStats["ISR"]["CALIBDIST"].items():
                for level in self.config.expectedDistributionLevels:
                    key = f"LSST CALIB {self.stageName.upper()} {ampName} DISTRIBUTION {level}-PCT"
                    row[ampName][f"darkDistribution_{level}"] = stats[key]

            # Append to output
            for ampName, stats in row.items():
                rowList.append(stats)

        return pipeBase.Struct(
            rowList=rowList,
        )


class CpVerifyRepackPhysicalFilterConnections(pipeBase.PipelineTaskConnections,
                                              dimensions={"instrument", "physical_filter"},
                                              defaultTemplate={}):
    """Connections class for calibration statistics with physical_filter
    (and instrument) dimensions.
    """
    detectorStats = cT.Input(
        name="detectorStats",
        doc="Input detector statistics.",
        storageClass="StructuredDataDict",
        dimensions={"instrument", "exposure", "detector"},
        multiple=True,
    )
    exposureStats = cT.Input(
        name="exposureStats",
        doc="Input exposure statistics.",
        storageClass="StructuredDataDict",
        dimensions={"instrument", "exposure"},
        multiple=True,
    )
    runStats = cT.Input(
        name="runStats",
        doc="Input Run statistics.",
        storageClass="StructuredDataDict",
        dimensions={"instrument"},
        multiple=True,
    )

    outputCatalog = cT.Output(
        name="cpvCatalog",
        doc="Output merged catalog.",
        storageClass="ArrowAstropy",
        dimensions={"instrument", "physical_filter"},
    )


class CpVerifyRepackPhysicalFilterConfig(pipeBase.PipelineTaskConfig,
                                         pipelineConnections=CpVerifyRepackPhysicalFilterConnections):
    expectedDistributionLevels = pexConfig.ListField(
        dtype=float,
        doc="Percentile levels expected in the calibration header.",
        default=[0, 5, 16, 50, 84, 95, 100],
    )


class CpVerifyRepackFlatTask(CpVerifyRepackTask):
    ConfigClass = CpVerifyRepackPhysicalFilterConfig

    stageName = "flat"

    def repackDetStats(self, detectorStats, detectorDims):
        rowList = []

        for detStats, detDims in zip(detectorStats, detectorDims):
            row = {}
            instrument = detDims["instrument"]
            exposure = detDims["exposure"]
            detector = detDims["detector"]

            # Get amp stats
            # AMP {ampName} [MEAN NOISE] value
            for ampName, stats in detStats["AMP"].items():
                row[ampName] = {
                    "instrument": instrument,
                    "exposure": exposure,
                    "detector": detector,
                    "amplifier": ampName,
                    "flatMean": stats["MEAN"],
                    "flatNoise": stats["NOISE"],
                }
            # Get catalog stats CATALOG

            # Get metadata stats
            # METADATA (RESIDUAL STDEV) {ampName} value

            # Get verify stats
            for ampName, stats in detStats["VERIFY"]["AMP"].items():
                row[ampName]["flatVerifyNoise"] = stats["NOISE"]
            # Get isr stats
            for ampName, stats in detStats["ISR"]["CALIBDIST"].items():
                for level in self.config.expectedDistributionLevels:
                    key = f"LSST CALIB {self.stageName.upper()} {ampName} DISTRIBUTION {level}-PCT"
                    row[ampName][f"flatDistribution_{level}"] = stats[key]
            # Get detector stats
            # DET
            row["detector"] = {"instrument": instrument,
                               "exposure": exposure,
                               "detector": detector,
                               "flatDetMean": detStats["DET"]["MEAN"],
                               "flatDetScatter": detStats["DET"]["SCATTER"],
                               }

            # Append to output
            for ampName, stats in row.items():
                rowList.append(stats)

        return pipeBase.Struct(
            rowList=rowList,
        )


class CpVerifyRepackDefectTask(CpVerifyRepackTask):
    stageName = "defects"

    def repackDetStats(self, detectorStats, detectorDims):
        rowList = []
        for detStats, detDims in zip(detectorStats, detectorDims):
            row = {}
            instrument = detDims["instrument"]
            detector = detDims["detector"]

            # Get amp stats
            for ampName, stats in detStats["AMP"].items():
                row[ampName] = {
                    "instrument": instrument,
                    "detector": detector,
                    "amplifier": ampName,
                }
            # Get catalog stats CATALOG
            # Get metadata stats METADATA
            # Get verify stats VERIFY
            # Get isr stats ISR
            nBadColumns = np.nan
            for ampName, stats in detStats["ISR"]["CALIBDIST"].items():
                if ampName == "detector":
                    nBadColumns = stats[ampName]["LSST CALIB DEFECTS N_BAD_COLUMNS"]
                else:
                    key = f"LSST CALIB DEFECTS {ampName} N_HOT"
                    row[ampName]["hotPixels"] = stats[ampName][key]
                    key = f"LSST CALIB DEFECTS {ampName} N_COLD"
                    row[ampName]["coldPixels"] = stats[ampName][key]
            # Get detector stats DET
            row["detector"] = {"instrument": instrument,
                               "detector": detector,
                               "nBadColumns": nBadColumns,
                               }
            for ampName, stats in row.items():
                rowList.append(stats)

        return rowList


class CpVerifyRepackCtiTask(CpVerifyRepackTask):
    stageName = "cti"
    pass


class CpVerifyRepackBfkTask(CpVerifyRepackTask):
    stageName = "bfk"
    pass


class CpVerifyRepackNoExpConnections(pipeBase.PipelineTaskConnections,
                                     dimensions={"instrument"},
                                     defaultTemplate={}):
    detectorStats = cT.Input(
        name="detectorStats",
        doc="Input detector statistics.",
        storageClass="StructuredDataDict",
        dimensions={"instrument", "detector"},
        multiple=True,
    )
    runStats = cT.Input(
        name="runStats",
        doc="Input Run statistics.",
        storageClass="StructuredDataDict",
        dimensions={"instrument"},
        multiple=True,
    )

    outputCatalog = cT.Output(
        name="cpvCatalog",
        doc="Output merged catalog.",
        storageClass="ArrowAstropy",
        dimensions={"instrument"},
    )


class CpVerifyRepackNoExpConfig(pipeBase.PipelineTaskConfig,
                                pipelineConnections=CpVerifyRepackNoExpConnections):

    expectedDistributionLevels = pexConfig.ListField(
        dtype=float,
        doc="Percentile levels expected in the calibration header.",
        default=[0, 5, 16, 50, 84, 95, 100],
    )
    hasMatrixCatalog = pexConfig.Field(
        dtype=bool,
        doc="Will a matrix catalog be created?",
        default=False,
    )


class CpVerifyRepackNoExpTask(CpVerifyRepackTask):
    """Repack cpVerify statistics for analysis_tools.

    Version for "verifyCalib" style results, which have no exposure
    dimension.
    """
    ConfigClass = CpVerifyRepackNoExpConfig
    _DefaultName = "cpVerifyRepack"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        inputs["detectorDims"] = [dict(exp.dataId.required) for exp in inputRefs.detectorStats]
        inputs["exposureDims"] = []
        inputs["exposureStats"] = []

        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)


class CpVerifyRepackPtcTask(CpVerifyRepackNoExpTask):
    stageName = "ptc"

    def repackDetStats(self, detectorStats, detectorDims):
        rowList = []

        for detStats, detDims in zip(detectorStats, detectorDims):
            row = {}

            instrument = detDims["instrument"]
            detector = detDims["detector"]

            # Get amp stats
            for ampName, stats in detStats["AMP"].items():
                row[ampName] = {
                    "instrument": instrument,
                    "detector": detector,
                    "amplifier": ampName,
                    "ampGain": stats["AMP_GAIN"],
                    "ampNoise": stats["AMP_NOISE"],
                    "ptcGain": stats["PTC_GAIN"],
                    "ptcNoise": stats["PTC_NOISE"],
                    "ptcTurnoff": stats["PTC_TURNOFF"],
                    "ptcFitType": stats["PTC_FIT_TYPE"],
                    "ptcBfeA00": stats["PTC_BFE_A00"],
                    "ptcRowMeanVariance": stats["PTC_ROW_MEAN_VARIANCE"],
                    "ptcRowMeanVarianceSlope": stats["PTC_ROW_MEAN_VARIANCE_SLOPE"],
                    "ptcMaxRawMeans": stats["PTC_MAX_RAW_MEANS"],
                    "ptcRawMeans": stats["PTC_RAW_MEANS"],
                    "ptcExpIdmask": stats["PTC_EXP_ID_MASK"],
                    "ptcCov10": stats["PTC_COV_10"],
                    "ptcCov10FitSlope": stats["PTC_COV_10_FIT_SLOPE"],
                    "ptcCov10FitOffset": stats["PTC_COV_10_FIT_OFFSET"],
                    "ptcCov10FitSuccess": stats["PTC_COV_10_FIT_SUCCESS"],
                    "ptcCov01": stats["PTC_COV_01"],
                    "ptcCov01FitSlope": stats["PTC_COV_01_FIT_SLOPE"],
                    "ptcCov01FitOffset": stats["PTC_COV_01_FIT_OFFSET"],
                    "ptcCov01FitSuccess": stats["PTC_COV_01_FIT_SUCCESS"],
                    "ptcCov11": stats["PTC_COV_11"],
                    "ptcCov11FitSlope": stats["PTC_COV_11_FIT_SLOPE"],
                    "ptcCov11FitOffset": stats["PTC_COV_11_FIT_OFFSET"],
                    "ptcCov11FitSuccess": stats["PTC_COV_11_FIT_SUCCESS"],
                    "ptcCov20": stats["PTC_COV_20"],
                    "ptcCov20FitSlope": stats["PTC_COV_20_FIT_SLOPE"],
                    "ptcCov20FitOffset": stats["PTC_COV_20_FIT_OFFSET"],
                    "ptcCov20FitSuccess": stats["PTC_COV_20_FIT_SUCCESS"],
                    "ptcCov02": stats["PTC_COV_02"],
                    "ptcCov02FitSlope": stats["PTC_COV_02_FIT_SLOPE"],
                    "ptcCov02FitOffset": stats["PTC_COV_02_FIT_OFFSET"],
                    "ptcCov02FitSuccess": stats["PTC_COV_02_FIT_SUCCESS"],
                }
            # Get catalog stats
            # Get detector stats
            # Get metadata stats
            # Get verify stats
            for ampName, stats in detStats["VERIFY"]["AMP"].items():
                row[ampName]["ptcVerifyGain"] = stats["PTC_GAIN"]
                row[ampName]["ptcVerifyNoise"] = stats["PTC_NOISE"]
                row[ampName]["ptcVerifyTurnoff"] = stats["PTC_TURNOFF"]
                row[ampName]["ptcVerifyBfeA00"] = stats["PTC_BFE_A00"]

            # Get isr stats

            # Append to output
            for ampName, stats in row.items():
                rowList.append(stats)

        return pipeBase.Struct(
            rowList=rowList,
        )


class CpVerifyRepackLinearityTask(CpVerifyRepackTask):
    stageName = "linearity"

    def repackDetStats(self, detectorStats, detectorDims):
        rowList = []

        for detStats, detDims in zip(detectorStats, detectorDims):
            row = {}

            instrument = detDims["instrument"]
            detector = detDims["detector"]

            # Get amp stats
            for ampName, stats in detStats["AMP"].items():
                centers, values = np.split(stats["LINEARITY_COEFFS"], 2)
                row[ampName] = {
                    "instrument": instrument,
                    "detector": detector,
                    "amplifier": ampName,
                    "fitParams": stats["FIT_PARAMS"],
                    "fitParamsErr": stats["FIT_PARAMS_ERR"],
                    "fitResiduals": stats["FIT_RESIDUALS"],
                    "splineCenters": centers,
                    "splineValues": values,
                    "linearityType": stats["LINEARITY_TYPE"],
                    "linearFit": stats["LINEAR_FIT"],
                }
            # Get catalog stats
            # Get detector stats
            # Get metadata stats
            # Get verify stats; no need to loop here.
            stats = detStats["VERIFY"]
            row["detector"] = {
                "instrument": instrument,
                "detector": detector,
                "linearityMaxResidualError": stats["MAX_RESIDUAL_ERROR"],
            }
            # Get isr stats

            # Append to output
            for ampName, stats in row.items():
                rowList.append(stats)

        return pipeBase.Struct(
            rowList=rowList,
        )


class CpVerifyRepackCrosstalkTask(CpVerifyRepackTask):
    pass
