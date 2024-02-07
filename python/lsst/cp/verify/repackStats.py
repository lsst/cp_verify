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

from astropy.table import Table

import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
import lsst.pex.config as pexConfig

__all__ = [
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
        dimensions={"instrument"},
        multiple=True,
    )

    outputCatalog = cT.Output(
        name="cpvCatalog",
        doc="Output merged catalog.",
        storageClass="ArrowAstropy",
        dimensions={"instrument"},
    )


class CpVerifyRepackInstrumentConfig(pipeBase.PipelineTaskConfig,
                                     pipelineConnections=CpVerifyRepackInstrumentConnections):

    expectedDistributionLevels = pexConfig.ListField(
        dtype=float,
        doc="Percentile levels expected in the calibration header.",
        default=[0, 5, 16, 50, 84, 95, 100],
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


class CpVerifyRepackDarkTask(CpVerifyRepackTask):
    stageName = "dark"


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


class CpVerifyRepackDefectTask(CpVerifyRepackTask):
    stageName = "defects"


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

        inputs["detectorDims"] = [dict(exp.dataId.required) for exp in inputRefs.detectorStats]
        inputs["exposureDims"] = []
        inputs["exposureStats"] = []

        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)


class CpVerifyRepackPtcTask(CpVerifyRepackNoExpTask):
    stageName = "ptc"


class CpVerifyRepackLinearityTask(CpVerifyRepackTask):
    stageName = "linearity"


class CpVerifyRepackCrosstalkTask(CpVerifyRepackTask):
    pass
