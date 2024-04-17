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
from astropy.table import vstack, Column, Table

import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
import lsst.pex.config as pexConfig


__all__ = ['CpVerifyExpMergeConfig', 'CpVerifyExpMergeTask',
           'CpVerifyRunMergeConfig', 'CpVerifyRunMergeTask',
           'CpVerifyExpMergeByFilterConfig', 'CpVerifyExpMergeByFilterTask',
           'CpVerifyRunMergeByFilterConfig', 'CpVerifyRunMergeByFilterTask',
           'CpVerifyVisitExpMergeConfig', 'CpVerifyVisitExpMergeTask',
           'CpVerifyVisitRunMergeConfig', 'CpVerifyVisitRunMergeTask',
           'CpVerifyCalibMergeConfig', 'CpVerifyCalibMergeTask']


class CpVerifyExpMergeConnections(pipeBase.PipelineTaskConnections,
                                  dimensions={"instrument", "exposure"},
                                  defaultTemplates={}):
    inputStats = cT.Input(
        name="detectorStats",
        doc="Input statistics to merge.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "exposure", "detector"],
        multiple=True,
    )
    inputResults = cT.Input(
        name="detectorResults",
        doc="Input results to merge.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "exposure", "detector"],
        multiple=True,
    )
    inputMatrix = cT.Input(
        name="detectorMatrix",
        doc="Input matrix to merge.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "exposure", "detector"],
        multiple=True,
    )
    camera = cT.PrerequisiteInput(
        name="camera",
        storageClass="Camera",
        doc="Input camera.",
        dimensions=["instrument", ],
        isCalibration=True,
    )

    outputStats = cT.Output(
        name="exposureStats",
        doc="Output statistics.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "exposure"],
    )
    outputResults = cT.Output(
        name="exposureResults",
        doc="Output results.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "exposure"],
    )
    outputMatrix = cT.Output(
        name="exposureMatrix",
        doc="Output matrix.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "exposure"],
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not self.config.hasMatrixCatalog:
            self.inputs.remove("inputMatrix")
            self.outputs.remove("outputMatrix")
        if not self.config.hasInputResults:
            self.inputs.remove("inputResults")


class CpVerifyExpMergeConfig(pipeBase.PipelineTaskConfig,
                             pipelineConnections=CpVerifyExpMergeConnections):
    """Configuration parameters for exposure stats merging.
    """
    statKeywords = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Dictionary of statistics to run on the set of detector values. The key should be the test "
        "name to record in the output, and the value should be the `lsst.afw.math` statistic name string.",
        default={},
    )
    hasMatrixCatalog = pexConfig.Field(
        dtype=bool,
        doc="Is there matrix catalog to merge?",
        default=False,
    )
    hasInputResults = pexConfig.Field(
        dtype=bool,
        doc="Are there results tables to merge?",
        default=False,
    )

    mergeDimension = pexConfig.Field(
        dtype=str,
        doc="Dimension name that these inputs will be merged over.",
        default="detector",
    )
    stageName = pexConfig.Field(
        dtype=str,
        doc="Stage name to use in any further analysis.",
        default="stageName",
    )


class CpVerifyExpMergeTask(pipeBase.PipelineTask):
    """Merge statistics from detectors together.
    """
    ConfigClass = CpVerifyExpMergeConfig
    _DefaultName = 'cpVerifyExpMerge'

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        dimensions = [dict(exp.dataId.required) for exp in inputRefs.inputStats]
        inputs['inputDims'] = dimensions

        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def run(self, inputStats, inputDims, camera=None, inputResults=None, inputMatrix=None,):
        """Merge statistics.

        Parameters
        ----------
        inputStats : `dict`
            A nested dictionary of measured statistics and
            verification results.
        inputDims : `dict`
            The input dimensions for each element of the inputStats.
        camera : `lsst.afw.cameraGeom.Camera`, optional
            The camera definition, used for identifying amplifier and
            detector names.
        inputResults : `astropy.Table`, optional
            The statistics information, formatted into a flat table.
        inputMatrix : `astropy.Table`, optional
            A table of results that represent the elements of matrices
            of values.

        Returns
        -------
        outputStats : `dict`
            A nested dictionary of merged statistics and verification
            results.
        outputResults : `astropy.Table`
            Flat table containing the merged results from all
            inputResults.
        outputMatrix : `astropy.Table`
            Table containing the merge results from all inputMatrix.

        See Also
        --------
        lsst.cp.verify.CpVerifyStatsTask
        """
        outputStats = {}   # This contains failure information
        success = True

        mergedStats = {}   # This contains the merged set of subcomponent stats.
        for inStats, dimensions in zip(inputStats, inputDims):
            thisId = dimensions[self.config.mergeDimension]
            thisName = thisId

            if self.config.mergeDimension == 'detector':
                thisName = camera[thisId].getName()

            calcStats = {}

            mergedStats[thisName] = inStats

            if inStats['SUCCESS'] is True:
                calcStats['SUCCESS'] = True
            else:
                calcStats['SUCCESS'] = False
                calcStats['FAILURES'] = list()
                success = False

                # See if we have verify information to check:
                if 'VERIFY' in inStats:
                    # See if an exposure failed
                    if 'EXP' in inStats['VERIFY'] and len(inStats['VERIFY']['EXP']) > 0:
                        expSuccess = inStats['VERIFY']['EXP'].pop('SUCCESS', False)
                        if not expSuccess:
                            for testName, testResult in inStats['VERIFY']['EXP'].items():
                                if testResult is False:
                                    calcStats['FAILURES'].append(testName)

                    # See if a detector failed
                    if 'DET' in inStats['VERIFY'] and len(inStats['VERIFY']['DET']) > 0:
                        detSuccess = inStats['VERIFY']['DET'].pop('SUCCESS', False)
                        if not detSuccess:
                            for testName, testResult in inStats['VERIFY']['DET'].items():
                                if testResult is False:
                                    calcStats['FAILURES'].append(testName)

                    # See if an amplifier failed
                    if 'AMP' in inStats['VERIFY'] and len(inStats['VERIFY']['AMP']) > 0:
                        for ampName, ampStats in inStats['VERIFY']['AMP'].items():
                            ampSuccess = ampStats.pop('SUCCESS')
                            if not ampSuccess:
                                for testName, testResult in ampStats.items():
                                    if testResult is False:
                                        calcStats['FAILURES'].append(ampName + " " + testName)

                    # See if a catalog failed
                    if 'CATALOG' in inStats['VERIFY'] and len(inStats['VERIFY']['CATALOG']) > 0:
                        for testName, testResult in inStats['VERIFY']['CATALOG'].items():
                            if testResult is False:
                                calcStats['FAILURES'].append(testName)
                else:
                    # No VERIFY info?  This must be partially accumulated.
                    # But we know there are failures somewhere.
                    # Drop any "SUCCESS" keys
                    _ = inStats.pop("SUCCESS", False)
                    for statKey, statDict in inStats.items():
                        if 'SUCCESS' in statDict and not statDict['SUCCESS']:
                            for failure in statDict['FAILURES']:
                                calcStats['FAILURES'].append(f"{statKey} {failure}")

            outputStats[thisName] = calcStats

        if self.config.mergeDimension == 'detector':
            outKey = 'EXP'
        else:
            outKey = 'RUN'

        groupSuccess = True
        if len(self.config.statKeywords):
            outputStats[outKey] = self.calcStatistics(mergedStats)
            outputStats['VERIFY'], groupSuccess = self.verify(mergedStats, outputStats)

        outputStats['SUCCESS'] = success & groupSuccess

        additionalResults = None
        if outKey in outputStats:
            # This is the only new information generated here.
            additionalResults, _ = self.pack(outputStats, inputDims, outKey)

        outputResults = self.mergeTable(inputResults, additionalResults)
        if inputMatrix is not None:
            outputMatrix = self.mergeTable(inputMatrix)
        else:
            outputMatrix = None

        return pipeBase.Struct(
            outputStats=outputStats,
            outputResults=outputResults,
            outputMatrix=outputMatrix,
        )

    @staticmethod
    def mergeTable(inputResults, newStats=None):
        """Merge input tables, padding columns as needed.

        Parameters
        ----------
        inputResults : `list` [`astropy.table.Table`]
            Input tables to merge.
        newStats : `astropy.table.Table`
            Additional table to merge.

        Returns
        -------
        merged : `astropy.table.Table`
            "Outer-join" merged table.
        """
        if inputResults is None:
            return Table()

        testTable = inputResults[0]  # This has the default set of columns.
        defaults = {key: -1 for key in testTable.columns}

        # Identify vector columns:
        for column in defaults.keys():
            if len(testTable[column].shape) > 1:
                defaults[column] = max(defaults[column],
                                       *[table[column].shape[1] for table in inputResults])

        # Pad vectors shorter than this:
        for column, length in defaults.items():
            if length > -1:  # is a vector
                for table in inputResults:
                    tableLength = table[column].shape[1]
                    if tableLength < length:  # this table is short
                        newColumn = []
                        for row in table[column]:
                            newColumn.append(np.pad(row,
                                                    (0, length-tableLength),
                                                    constant_values=np.nan))
                            table[column] = Column(newColumn, name=column, unit=table[column].unit)

        if len(inputResults) > 0:
            outputResults = vstack(inputResults)
        else:
            outputResults = inputResults

        if newStats:
            return vstack([outputResults, Table(newStats)])
        else:
            return outputResults

    def calcStatistics(self, statisticsDict):
        """Calculate exposure level statistics based on the existing
        per-amplifier and per-detector measurements.

        Parameters
        ----------
        statisticsDictionary : `dict [`str`, `dict` [`str`, scalar]],
            Dictionary of measured statistics.  The top level
            dictionary is keyed on the detector names, and contains
            the measured statistics from the per-detector
            measurements.

        Returns
        -------
        outputStatistics : `dict` [`str, scalar]
            A dictionary of the statistics measured and their values.
        """
        raise NotImplementedError("Subclasses must implement verification criteria.")

    def verify(self, detectorStatistics, statisticsDictionary):

        """Verify if the measured statistics meet the verification criteria.

        Parameters
        ----------
        detectorStatistics : `dict` [`str`, `dict` [`str`, scalar]]
            Merged set of input detector level statistics.
        statisticsDictionary : `dict` [`str`, `dict` [`str`, scalar]]
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

        Raises
        ------
        NotImplementedError :
            This method must be implemented by the calibration-type
            subclass.
        """
        raise NotImplementedError("Subclasses must implement verification criteria.")


class CpVerifyRunMergeConnections(pipeBase.PipelineTaskConnections,
                                  dimensions={"instrument", },
                                  defaultTemplates={}):
    inputStats = cT.Input(
        name="exposureStats",
        doc="Input statistics to merge.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "exposure"],
        multiple=True,
    )
    inputResults = cT.Input(
        name="exposureResults",
        doc="Input results table to merge.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "exposure"],
        multiple=True,
    )
    inputMatrix = cT.Input(
        name="exposureMatrix",
        doc="Input matrix table to merge.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "exposure"],
        multiple=True,
    )
    camera = cT.PrerequisiteInput(
        name="camera",
        storageClass="Camera",
        doc="Input camera.",
        dimensions=["instrument", ],
        isCalibration=True,
    )

    outputStats = cT.Output(
        name="runStats",
        doc="Output statistics.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", ],
    )
    outputResults = cT.Output(
        name="runResults",
        doc="Output merged results table.",
        storageClass="ArrowAstropy",
        dimensions=["instrument",],
    )
    outputMatrix = cT.Output(
        name="runMatrix",
        doc="Output merged matrix table.",
        storageClass="ArrowAstropy",
        dimensions=["instrument",],
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not self.config.hasMatrixCatalog:
            self.inputs.remove("inputMatrix")
            self.outputs.remove("outputMatrix")


class CpVerifyRunMergeConfig(CpVerifyExpMergeConfig,
                             pipelineConnections=CpVerifyRunMergeConnections):
    """Configuration paramters for exposure stats merging.
    """
    mergeDimension = pexConfig.Field(
        dtype=str,
        doc="Dimension name for this input.",
        default="exposure",
    )


class CpVerifyRunMergeTask(CpVerifyExpMergeTask):
    """Merge statistics from detectors together.
    """
    ConfigClass = CpVerifyRunMergeConfig
    _DefaultName = 'cpVerifyRunMerge'

    pass
# End ExpMerge/RunMerge


# Begin ExpMergeByFilter/RunMergeByFilter
class CpVerifyExpMergeByFilterConnections(pipeBase.PipelineTaskConnections,
                                          dimensions={"instrument", "exposure", "physical_filter"},
                                          defaultTemplates={}):
    inputStats = cT.Input(
        name="exposureStats",
        doc="Input statistics to merge.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "exposure", "detector", "physical_filter"],
        multiple=True,
    )
    inputResults = cT.Input(
        name="exposureResults",
        doc="Input results table to merge.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "exposure", "detector", "physical_filter"],
        multiple=True,
    )
    inputMatrix = cT.Input(
        name="exposureMatrix",
        doc="Input matrix table to merge.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "exposure", "detector", "physical_filter"],
        multiple=True,
    )
    camera = cT.PrerequisiteInput(
        name="camera",
        storageClass="Camera",
        doc="Input camera.",
        dimensions=["instrument", ],
        isCalibration=True,
    )

    outputStats = cT.Output(
        name="runStats",
        doc="Output statistics.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "exposure", "physical_filter"],
    )
    outputResults = cT.Output(
        name="runResults",
        doc="Output merged results table.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "exposure", "physical_filter"],
    )
    outputMatrix = cT.Output(
        name="runMatrix",
        doc="Output merged matrix table.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "exposure", "physical_filter"],
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not self.config.hasMatrixCatalog:
            self.inputs.remove("inputMatrix")
            self.outputs.remove("outputMatrix")


class CpVerifyExpMergeByFilterConfig(CpVerifyExpMergeConfig,
                                     pipelineConnections=CpVerifyExpMergeByFilterConnections):
    """Configuration paramters for exposure stats merging.
    """
    mergeDimension = pexConfig.Field(
        dtype=str,
        doc="Dimension name for this input.",
        default="detector",
    )


class CpVerifyExpMergeByFilterTask(CpVerifyExpMergeTask):
    """Merge statistics from detectors together.
    """
    ConfigClass = CpVerifyRunMergeConfig
    _DefaultName = 'cpVerifyRunMerge'

    pass


class CpVerifyRunMergeByFilterConnections(pipeBase.PipelineTaskConnections,
                                          dimensions={"instrument", "physical_filter"},
                                          defaultTemplates={}):
    inputStats = cT.Input(
        name="exposureStats",
        doc="Input statistics to merge.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "exposure", "physical_filter"],
        multiple=True,
    )
    inputResults = cT.Input(
        name="exposureResults",
        doc="Input results table to merge.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "exposure", "physical_filter"],
        multiple=True,
    )
    inputMatrix = cT.Input(
        name="exposureMatrix",
        doc="Input matrix table to merge.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "exposure", "physical_filter"],
        multiple=True,
    )
    camera = cT.PrerequisiteInput(
        name="camera",
        storageClass="Camera",
        doc="Input camera.",
        dimensions=["instrument", ],
        isCalibration=True,
    )

    outputStats = cT.Output(
        name="runStats",
        doc="Output statistics.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "physical_filter"],
    )
    outputResults = cT.Output(
        name="runResults",
        doc="Output merged results table.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "physical_filter"],
    )
    outputMatrix = cT.Output(
        name="runMatrix",
        doc="Output merged matrix table.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "physical_filter"],
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not self.config.hasMatrixCatalog:
            self.inputs.remove("inputMatrix")
            self.outputs.remove("outputMatrix")


class CpVerifyRunMergeByFilterConfig(CpVerifyRunMergeConfig,
                                     pipelineConnections=CpVerifyRunMergeByFilterConnections):
    """Configuration paramters for exposure stats merging.
    """
    pass


class CpVerifyRunMergeByFilterTask(CpVerifyExpMergeTask):
    """Merge statistics from detectors together.
    """
    ConfigClass = CpVerifyRunMergeByFilterConfig
    _DefaultName = 'cpVerifyRunMergeByFilter'

    pass
# End ExpMergeByFilter/RunMergeByFilter


# Begin ExpMergeByVisit/RunMergeByVisit
class CpVerifyVisitExpMergeConnections(pipeBase.PipelineTaskConnections,
                                       dimensions={"instrument", "visit"},
                                       defaultTemplates={}):
    inputStats = cT.Input(
        name="detectorStats",
        doc="Input statistics to merge.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "visit", "detector"],
        multiple=True,
    )
    inputResults = cT.Input(
        name="detectorResults",
        doc="Input results to merge.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "visit", "detector"],
        multiple=True,
    )
    inputMatrix = cT.Input(
        name="detectorMatrix",
        doc="Input matrix to merge.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "visit", "detector"],
        multiple=True,
    )
    camera = cT.PrerequisiteInput(
        name="camera",
        storageClass="Camera",
        doc="Input camera.",
        dimensions=["instrument", ],
        isCalibration=True,
    )

    outputStats = cT.Output(
        name="exposureStats",
        doc="Output statistics.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "visit"],
    )
    outputResults = cT.Output(
        name="exposureResults",
        doc="Output results.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "visit"],
    )
    outputMatrix = cT.Output(
        name="exposureMatrix",
        doc="Output matrix.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "visit"],
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not self.config.hasMatrixCatalog:
            self.inputs.remove("inputMatrix")
            self.outputs.remove("outputMatrix")


class CpVerifyVisitExpMergeConfig(CpVerifyExpMergeConfig,
                                  pipelineConnections=CpVerifyVisitExpMergeConnections):
    pass


class CpVerifyVisitExpMergeTask(CpVerifyExpMergeTask):
    """Merge visit based data."""

    ConfigClass = CpVerifyVisitExpMergeConfig
    _DefaultName = 'cpVerifyVisitExpMerge'

    pass


class CpVerifyVisitRunMergeConnections(pipeBase.PipelineTaskConnections,
                                       dimensions={"instrument"},
                                       defaultTemplates={}):
    inputStats = cT.Input(
        name="exposureStats",
        doc="Input statistics to merge.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "visit"],
        multiple=True,
    )
    inputResults = cT.Input(
        name="exposureResults",
        doc="Input results table to merge.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "visit"],
        multiple=True,
    )
    inputMatrix = cT.Input(
        name="exposureMatrix",
        doc="Input matrix table to merge.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "visit"],
        multiple=True,
    )
    camera = cT.PrerequisiteInput(
        name="camera",
        storageClass="Camera",
        doc="Input camera.",
        dimensions=["instrument", ],
        isCalibration=True,
    )

    outputStats = cT.Output(
        name="runStats",
        doc="Output statistics.",
        storageClass="StructuredDataDict",
        dimensions=["instrument"],
    )
    outputResults = cT.Output(
        name="runResults",
        doc="Output merged results table.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", ],
    )
    outputMatrix = cT.Output(
        name="runMatrix",
        doc="Output merged matrix table.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", ],
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not self.config.hasMatrixCatalog:
            self.inputs.remove("inputMatrix")
            self.outputs.remove("outputMatrix")


class CpVerifyVisitRunMergeConfig(CpVerifyRunMergeConfig,
                                  pipelineConnections=CpVerifyVisitRunMergeConnections):
    pass


class CpVerifyVisitRunMergeTask(CpVerifyRunMergeTask):
    """Merge visit based data."""

    ConfigClass = CpVerifyVisitRunMergeConfig
    _DefaultName = 'cpVerifyVisitRunMerge'

    pass
# End ExpMergeByVisit/RunMergeByVisit


# Begin CalibMerge (this is a one-step)
class CpVerifyCalibMergeConnections(pipeBase.PipelineTaskConnections,
                                    dimensions={"instrument"},
                                    defaultTemplates={}):
    inputStats = cT.Input(
        name="exposureStats",
        doc="Input statistics to merge.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "detector"],
        multiple=True,
    )
    inputResults = cT.Input(
        name="exposureResults",
        doc="Input results table to merge.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "detector"],
        multiple=True,
    )
    inputMatrix = cT.Input(
        name="exposureMatrix",
        doc="Input matrix table to merge.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "detector"],
        multiple=True,
    )
    camera = cT.PrerequisiteInput(
        name="camera",
        storageClass="Camera",
        doc="Input camera.",
        dimensions=["instrument", ],
        isCalibration=True,
    )

    outputStats = cT.Output(
        name="exposureStats",
        doc="Output statistics.",
        storageClass="StructuredDataDict",
        dimensions=["instrument"],
    )
    outputResults = cT.Output(
        name="runResults",
        doc="Output merged results table.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", ],
    )
    outputMatrix = cT.Output(
        name="runMatrix",
        doc="Output merged matrix table.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", ],
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not self.config.hasMatrixCatalog:
            self.inputs.remove("inputMatrix")
            self.outputs.remove("outputMatrix")


class CpVerifyCalibMergeConfig(CpVerifyRunMergeConfig,
                               pipelineConnections=CpVerifyCalibMergeConnections):
    """Configuration paramters for exposure stats merging.
    """
    pass


class CpVerifyCalibMergeTask(CpVerifyRunMergeTask):
    """Merge statistics from detectors together.
    """
    ConfigClass = CpVerifyCalibMergeConfig
    _DefaultName = 'cpVerifyCalibMerge'

    pass
# End CalibMerge
