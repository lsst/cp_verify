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
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
import lsst.afw.image as afwImage
import pandas as pd
from lsst.analysis.tools.actions.plot import FocalPlaneGeometryPlot as FPGPlot
from .verifyDefects import (
    CpVerifyDefectsConnections,
)
from .verifyStats import (
    CpVerifyStatsConfig,
)

__all__ = ["CpMeasureDefectsStabilityTaskConfig", "CpMeasureDefectsStabilityTask"]


def lookupStaticCalibrations(datasetType, registry, quantumDataId, collections):
    # This will retrieve all the filter transmissions available
    # associated with the {instrument, detector} dimensions.
    timespan = Timespan(begin=None, end=None)
    result = []
    # First iterate over all of the data IDs for this dataset type that are
    # consistent with the quantum data ID.
    for dataId in registry.queryDataIds(datasetType.dimensions, dataId=quantumDataId):
        # Find the dataset with this data ID using the unbounded timespan.
        if ref := registry.findDataset(
            datasetType, dataId, collections=collections, timespan=timespan
        ):
            result.append(ref)
    return result


class CpMeasureDefectsStabilityConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("instrument", "detector"),  # TODO: Update the inheritance
):
    referenceDefects = cT.PrerequisiteInput(
        name="defects",
        doc="Reference defects.",
        storageClass="Defects",
        dimensions=(
            "instrument",
            "detector",
        ),
        lookupFunction=lookupStaticCalibrations,
        deferLoad=True,
        multiple=True,
    )

    productionDefects = cT.PrerequisiteInput(
        name="defects",
        doc="Production defects.",
        storageClass="Defects",
        dimensions=(
            "instrument",
            "detector",
        ),
        lookupFunction=lookupStaticCalibrations,
        deferLoad=True,
        multiple=True,
    )

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera associated with these defects.",
        storageClass="Camera",
        dimensions=("instrument",),
    )
    outputTable = cT.Output(
        name="defectsStability",
        doc="Dataframe storing the variation of defects.",
        storageClass="DataFrame",
        dimensions=(),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)


class CpMeasureDefectsStabilityTaskConfig(
    CpVerifyStatsConfig, pipelineConnections=CpVerifyDefectsConnections
):
    """
    Configuration for measuring stability
    of defects from individual defect masks.
    """

    refLabel = pexConfig.Field(
        dtype=str,
        doc="Label to apply to the reference defect mask.",
        optional=False,
    )

    prodLabel = pexConfig.Field(
        dtype=str,
        doc="Label to apply to the production defect mask.",
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
        self.detName_Obj = self.camera.getNameMap()
        self.det_amp = []
        for detectorName, detectorObject in self.detName_Obj.items():
            for amp in detectorObject.getAmplifiers():
                self.det_amp.append(f"{detectorName}_{amp.getName()}")

        """
        Here, setup the default dataframe based on the
        det_amp object and the camera object
        """
        columns = []
        for plane in self.referenceDefects[0].getMaskPlaneDict().keys():
            if plane not in self.skipMasks:
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


class CpMeasureDefectsStabilityTask(pipeBase.PipelineTask):
    """Task to measure variation  defects in combined images."""

    ConfigClass = CpMeasureDefectsStabilityTaskConfig
    _DefaultName = "cpMeasureDefectsStability"

    def evaluate_mask_variability(self, masks, dicts):
        """
        Evaluate the variability of each mask plane across an entire
        detector/amplifier

        masks: List of 2D numpy arrays, each corresponding to the
               reference and production defect mask for a given
               detector/amplifier.
        dicts: List of bit_position dictionaries corresponding to
               each mask.
        """
        # Find the union of all defect names across all masks
        all_defect_names = set().union(*(d.keys() for d in dicts))

        stacks = {}

        for name in all_defect_names:
            """
            Extract the specific binary plane for
            this defect from EVERY mask
            If a mask doesn't have this defect, it's treated
            as an empty plane (zeros)"""
            planes = []
            for i in range(len(masks)):
                m = masks[i]
                d = dicts[i]
                if name in d:
                    bit_pos = d[name]
                    binary_plane = (m & (1 << bit_pos)) > 0
                else:
                    binary_plane = np.zeros_like(m, dtype=bool)

                planes.append(binary_plane.astype(float))

            # Convert list to a 3D stack (N, Height, Width)
            stack = np.stack(planes)

            stacks[name] = stack
        return stacks

    def run(self, referenceDefects, productionDefects, camera):

        # Perform the task here, measuring the variance in the defects
        defectsStability = self.processDetectorsAndAmps(
            referenceDefects, productionDefects, self.config.resultDataFrame
        )

        return pipeBase.Struct(
            defectsStability=defectsStability,
        )

    def processDetectorsAndAmps(self, referenceDefects, productionDefects, defectData):
        for detectorName, detectorObject in self.config.detName_Obj.items():
            # Select identical detectors
            defect1 = referenceDefects[
                [
                    x.getMetadata()["DETECTOR"] == detectorObject.getId()
                    for x in referenceDefects
                ]
            ]
            defect2 = productionDefects[
                [
                    x.getMetadata()["DETECTOR"] == detectorObject.getId()
                    for x in productionDefects
                ]
            ]

            # Make two mask objects that possess the defect masks
            exp1 = afwImage.ExposureF(detectorObject.getBBox())
            defect1.maskPixels(exp1.maskedImage)
            mask1 = exp1.maskedImage.mask

            exp2 = afwImage.ExposureF(detectorObject.getBBox())
            defect2.maskPixels(exp2.maskedImage)
            mask2 = exp2.maskedImage.mask

            del exp1, exp2

            """
            # Quick check on the nature of the two defect masks
            # If they are identical, there's probably an
            # elegant way to speed this up
            if all((mask1.array - mask2.array).flatten() == 0):
                print("Mask planes are identical.")
            else:
                print("Mask planes are not identical.")
            # NOTE: Commenting the above out for now, there is definitely
            # a faster way to address identical mask planes...
            """

            # Assert that the masks have the same plane definitions
            assert (
                mask1.getMaskPlaneDict() == mask2.getMaskPlaneDict()
            ), "Mask plane dictionaries are incompatible"

            # For each amp
            for amp in detectorObject.getAmplifiers():

                if self.config.includeAny:
                    # Instantiate the anymask object
                    anyMasks = []
                    for acq in [self.config.refLabel, self.config.prodLabel]:
                        anyMasks.append(np.zeros(np.shape(mask1[amp.getBBox()].array)))

                # Get the amplifier name and compute the
                # variability inside that amplifier region
                ampName = amp.getName()
                variability = self.evaluate_mask_variability(
                    [mask1[amp.getBBox()].array, mask2[amp.getBBox()].array],
                    [mask1.getMaskPlaneDict(), mask2.getMaskPlaneDict()],
                )

                # For each mask plane
                for plane, planeVals in variability.items():
                    if plane not in self.config.skipMasks:
                        # Compute NPix for each run and add it to the df
                        for i, acq in enumerate(
                            [self.config.refLabel, self.config.prodLabel]
                        ):
                            defectData.loc[
                                f"{detectorName}_{ampName}", f"{acq}_{plane}"
                            ] = np.count_nonzero(planeVals[i])
                            # N pix that are this mask plane type, in
                            # this acq run, indexed by amplifier
                            anyMasks[i] += planeVals[i]

                        # Compute the difference and add it to the df
                        defectData.loc[
                            f"{detectorName}_{ampName}",
                            f"{self.config.acqString}{plane}",
                        ] = np.sum(planeVals[i - 1] - planeVals[i], dtype=int)
                        # Difference in N pix that are this
                        # mask plane type, indexed by amplifier

                if self.config.includeAny:
                    # Perform the same calculation,
                    # here for the anyMask case
                    for i, acq in enumerate(
                        [self.config.refLabel, self.config.prodLabel]
                    ):
                        defectData.loc[f"{detectorName}_{ampName}", f"{acq}_Any"] = (
                            np.count_nonzero(anyMasks[i])
                        )
                    defectData.loc[
                        f"{detectorName}_{ampName}", f"{self.config.acqString}Any"
                    ] = np.sum(anyMasks[i - 1] - anyMasks[i], dtype=int)

        return defectData

    # def runQuantum(
    #     self, butlerQC, inputRefs, outputRefs
    # ):  # TODO: rewrite this function, copied from MergeDefectsCombinedTask
    #     inputs = butlerQC.get(inputRefs)
    #     # Turn inputFlatDefects and inputDarkDefects into a list which
    #     # is what MergeDefectsTask expects.  If there are multiple,
    #     # use the one with the most inputs.
    #     tempList = [
    #         self.chooseBest(inputs["inputFlatDefects"]),
    #         self.chooseBest(inputs["inputDarkDefects"]),
    #         self.chooseBest(inputs["inputBiasDefects"]),
    #     ]

    #     if "inputManualDefects" in inputs.keys():
    #         tempList.extend(inputs["inputManualDefects"])

    #     # Rename inputDefects
    #     inputsCombined = {"inputDefects": tempList,
    #                       "camera": inputs["camera"]}

    #     outputs = super().run(**inputsCombined)
    #     butlerQC.put(outputs, outputRefs)


class CpPlotDefectsStabilityTaskConnections(
    pipeBase.PipelineTaskConnections, dimensions=()
):
    """
    Connections for plotting stability
    of defects from individual defect masks.
    """

    inputTable = cT.PrerequisiteInput(
        name="defectsStability",
        doc="Table of defect stability measurements.",
        storageClass="DataFrame",
        dimensions=(),
    )

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Camera associated with these defects.",
        storageClass="Camera",
        dimensions=("instrument",),
    )

    outputTable = cT.Output(
        name="defectStabilityPlots",
        doc="Plots showing the variation of defect planes.",
        storageClass="Plot",
        dimensions=("instrument",),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)


class CpPlotDefectsStabilityTaskConfig(
    CpVerifyStatsConfig
):  # TODO: Update this inheritance?
    """
    Configuration for plotting stability
    of defects from individual defect masks.
    """

    skipPlanes = pexConfig.ListField(
        dtype=str,
        doc="Mask planes to skip when plotting.",
        default=(),
        optional=True,
    )

    _filter = pexConfig.Field(
        dtype=str,
        doc="Filter that the defects were generated from, used in plotting.",
        default="None",
        optional=True,
    )

    tract = pexConfig.Field(
        dtype=str,
        doc="Tract that the defects were generated from, used in plotting.",
        default="None",
        optional=True,
    )

    skymap = pexConfig.Field(
        dtype=str,
        doc="Skymap that the defects were generated from, used in plotting.",
        default="None",
        optional=True,
    )

    bands = pexConfig.ListField(
        dtype=str,
        doc="Bands that the defects were generated from, used in plotting.",
        default="None",
        optional=True,
    )

    plotName = pexConfig.Field(
        dtype=str,
        doc="Plot name.",
        default=None,
        optional=True,
    )

    tableName = pexConfig.Field(
        dtype=str,
        doc="Table name.",
        default=None,
        optional=True,
    )

    plotMax = pexConfig.Field(
        dtype=float,
        doc="Maximum value to be shown in the plot.",
        default=500,
        optional=True,
    )

    plotMin = pexConfig.Field(
        dtype=float,
        doc="Minimum value to be shown in the plot.",
        default=-500,
        optional=True,
    )

    def validate(self):
        super().validate()


class CpPlotDefectsStabilityTask(pipeBase.PipelineTask):
    """
    Task for plotting stability of defects
    from individual defect masks.
    """

    ConfigClass = CpPlotDefectsStabilityTaskConfig
    _DefaultName = "CpPlotDefectsStability"

    def prepDF(self, inputData, col):
        """
        Prepare the dataframe to be used by the FocalPlaneGeometryPlot
        """

        # df is indexed by RXX_SYY_CZZ
        amps = [x[8:] for x in inputData.index]
        dets = [x[:7] for x in inputData.index]

        columns = ["detector", "amplifier", "z"]
        data = np.array([dets, amps, inputData[col].tolist()]).T
        retDF = pd.DataFrame(data=data, columns=columns, dtype=str)
        retDF["z"] = retDF["z"].astype("float32")

        return retDF

    def run(self, inputTable):

        figs = self.makePlots(inputTable)

        return pipeBase.Struct(
            defectStabilityPlots=figs,
        )

    def makePlots(self, inputTable):
        figs = []
        fpgObj = FPGPlot()
        plotInfo = {
            "run": self.config.run,
            "skymap": self.config.run,
            "filter": self.config._filter,
            "tract": self.config.tract,
            "bands": self.config.bands,
            "plotName": self.config.plotName,
            "tableName": self.config.tableName,
        }
        for col in inputTable.columns.values:
            if col not in self.config.skipPlanes:
                for key, cfig in zip(
                    ["plotName", "tableName"],
                    [self.config.plotName, self.config.tableName],
                ):
                    if cfig is None:
                        plotInfo[key] = f"{col}"
                subTable = self.prepDF(inputTable, col)
                figs.append(
                    fpgObj.makePlot(
                        subTable,
                        self.camera,
                        plotInfo,
                        plotMax=self.config.plotMax,
                        plotMin=self.config.plotMin,
                    )
                )
        return figs
