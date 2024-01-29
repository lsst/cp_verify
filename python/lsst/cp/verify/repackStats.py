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
    "CpVerifyRepackConnections",
    "CpVerifyRepackConfig",
    "CpVerifyRepackBiasTask",
    "CpVerifyRepackDarkTask",
    "CpVerifyRepackFlatTask",
    "CpVerifyRepackDefectTask",
    "CpVerifyRepackPtcTask",
    "CpVerifyRepackBfkTask",
    "CpVerifyRepackCtiTask",
]


class CpVerifyRepackConnections(pipeBase.PipelineTaskConnections,
                                dimensions={"instrument"},
                                defaultTemplate={}):
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
        dimensions={"instrument"},
    )


class CpVerifyRepackConfig(pipeBase.PipelineTaskConfig,
                           pipelineConnections=CpVerifyRepackConnections):

    expectedDistributionLevels = pexConfig.ListField(
        dtype=float,
        doc="Percentile levels expected in the calibration header.",
        default=[0, 5, 16, 50, 84, 95, 100],
    )


class CpVerifyRepackTask(pipeBase.PipelineTask):
    """Repack cpVerify statistics for analysis_tools.
    """
    ConfigClass = CpVerifyRepackConfig
    _DefaultName = "cpVerifyRepack"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        inputs["detectorDims"] = [exp.dataId.byName() for exp in inputRefs.detectorStats]
        inputs["exposureDims"] = [exp.dataId.byName() for exp in inputRefs.exposureStats]

        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def run(self, detectorStats, detectorDims, exposureStats, exposureDims, runStats):
        """
        """
        rowList = self.repack(detectorStats, detectorDims, exposureStats, exposureDims, runStats)
        catalog = Table(rowList)

        return pipeBase.Struct(
            outputCatalog=catalog,
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
                row[ampName]["biasSerialProfile"] = projStats["AMP_HPROJECTION"][ampName]
            for ampName in projStats["AMP_VPROJECTION"].keys():
                row[ampName]["biasParallelProfile"] = projStats["AMP_VPROJECTION"][ampName]

            # Create output table:
            for ampName, stats in row.items():
                rowList.append(stats)
        return rowList


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

        return rowList


class CpVerifyRepackFlatTask(CpVerifyRepackTask):
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

        return rowList


class CpVerifyRepackDefectTask(CpVerifyRepackTask):
    stageName = "defects"
    pass


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


class CpVerifyRepackNoExpTask(CpVerifyRepackTask):
    """Repack cpVerify statistics for analysis_tools.

    Version for "verifyCalib" style results, which have no exposure
    dimension.
    """
    ConfigClass = CpVerifyRepackNoExpConfig
    _DefaultName = "cpVerifyRepack"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        inputs["detectorDims"] = [exp.dataId.byName() for exp in inputRefs.detectorStats]
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
                    "ptcRawMeans": stats["PTC_RAW_MEANS"],
                    "ptcExpIdmask": stats["PTC_EXP_ID_MASK"],
                    "ptcCov10": stats["PTC_COV_10"],
                    "ptcCov01": stats["PTC_COV_01"],
                    "ptcCov11": stats["PTC_COV_11"],
                    "ptcCov20": stats["PTC_COV_20"],
                    "ptcCov02": stats["PTC_COV_02"],
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

        return rowList


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

        return rowList


class CpVerifyRepackCrosstalkTask(CpVerifyRepackTask):
    pass
