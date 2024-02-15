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
from astropy.table import vstack, Column

import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
import lsst.pex.config as pexConfig


__all__ = ['CpVerifyExpMergeConfig', 'CpVerifyExpMergeTask',
           'CpVerifyRunMergeConfig', 'CpVerifyRunMergeTask',
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


class CpVerifyExpMergeConfig(pipeBase.PipelineTaskConfig,
                             pipelineConnections=CpVerifyExpMergeConnections):
    """Configuration parameters for exposure stats merging.
    """
    exposureStatKeywords = pexConfig.DictField(
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

    def run(self, inputStats, camera, inputDims=None, inputResults=None, inputMatrix=None,):
        """Merge statistics.

        Parameters
        ----------
        inputStats : `list` [`dict`]
            Measured statistics for a detector (from
            CpVerifyStatsTask).
        camera : `lsst.afw.cameraGeom.Camera`
            The camera geometry for this exposure.
        inputDims : `list` [`dict`]
            List of dictionaries of input data dimensions/values.
            Each list entry should contain:

            ``"exposure"``
                exposure id value (`int`)
            ``"detector"``
                detector id value (`int`)

        Returns
        -------
        outputStats

        See Also
        --------
        lsst.cp.verify.CpVerifyStatsTask

        Notes
        -----
        The outputStats should have a yaml representation of the form:

        DET:
          DetName1:
            FAILURES:
              - TEST_NAME
            STAT: value
            STAT2: value2
          DetName2:
        VERIFY:
           TEST: boolean
           TEST2: boolean
        SUCCESS: boolean
        """
        outputStats = {}
        success = True

        mergedStats = {}
        for detStats, dimensions in zip(inputStats, inputDims):
            detId = dimensions['detector']
            detName = camera[detId].getName()
            calcStats = {}

            mergedStats[detName] = detStats

            if detStats['SUCCESS'] is True:
                calcStats['SUCCESS'] = True
            else:
                calcStats['SUCCESS'] = False
                calcStats['FAILURES'] = list()
                success = False
                # See if the detector failed
                if 'DET' in detStats['VERIFY']:
                    detSuccess = detStats['VERIFY']['DET'].pop('SUCCESS', False)
                    if not detSuccess:
                        for testName, testResult in detStats['VERIFY']['DET'].items():
                            if testResult is False:
                                calcStats['FAILURES'].append(testName)
                # See if the catalog failed
                if 'CATALOG' in detStats['VERIFY']:
                    for testName, testResult in detStats['VERIFY']['CATALOG'].items():
                        if testResult is False:
                            calcStats['FAILURES'].append(testName)
                # See if an amplifier failed
                for ampName, ampStats in detStats['VERIFY']['AMP'].items():
                    ampSuccess = ampStats.pop('SUCCESS')
                    if not ampSuccess:
                        for testName, testResult in ampStats.items():
                            if testResult is False:
                                calcStats['FAILURES'].append(ampName + " " + testName)

            outputStats[detName] = calcStats

        exposureSuccess = True
        if len(self.config.exposureStatKeywords):
            outputStats['EXP'] = self.exposureStatistics(mergedStats)
            outputStats['VERIFY'], exposureSuccess = self.verify(mergedStats, outputStats)

        outputStats['SUCCESS'] = success & exposureSuccess

        outputResults = mergeTable(inputResults, outputStats)
        outputMatrix = mergeTable(inputMatrix, outputStats)
        return pipeBase.Struct(
            outputStats=outputStats,
            outputResults=outputResults,
            outputMatrix=outputMatrix,
        )

    def exposureStatistics(self, statisticsDict):
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
                                  dimensions={"instrument"},
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


class CpVerifyRunMergeConfig(pipeBase.PipelineTaskConfig,
                             pipelineConnections=CpVerifyRunMergeConnections):
    """Configuration paramters for exposure stats merging.
    """
    runStatKeywords = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Dictionary of statistics to run on the set of exposure values. The key should be the test "
        "name to record in the output, and the value should be the `lsst.afw.math` statistic name string.",
        default={},
    )
    hasMatrixCatalog = pexConfig.Field(
        dtype=bool,
        doc="Is there matrix catalog to merge?",
        default=False,
    )


class CpVerifyRunMergeTask(pipeBase.PipelineTask):
    """Merge statistics from detectors together.
    """
    ConfigClass = CpVerifyRunMergeConfig
    _DefaultName = 'cpVerifyRunMerge'

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        dimensions = [dict(exp.dataId.required) for exp in inputRefs.inputStats]
        inputs['inputDims'] = dimensions

        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def run(self, inputStats, inputDims, inputResults=None, inputMatrix=None):
        """Merge statistics.

        Parameters
        ----------
        inputStats : `list` [`dict`]
            Measured statistics for a detector.
        inputDims : `list` [`dict`]
            List of dictionaries of input data dimensions/values.
            Each list entry should contain:

            ``"exposure"``
                exposure id value (`int`)

        inputResults : `list` [`astropy.table.Table`]
            List of input tables of results to merge.
        inputMatrix : `list` [`astropy.table.Table`]
            List of input matrix tables to merge.

        Returns
        -------
        outputStats : 


        Notes
        -----
        The outputStats should have a yaml representation as follows.

        VERIFY:
          ExposureId1:
            VERIFY_TEST1: boolean
            VERIFY_TEST2: boolean
          ExposureId2:
            [...]
          TEST_VALUE: boolean
          TEST_VALUE2: boolean
        """
        outputStats = {}
        success = True
        for expStats, dimensions in zip(inputStats, inputDims):
            expId = dimensions.get('exposure', dimensions.get('visit', None))
            if expId is None:
                raise RuntimeError("Could not identify the exposure from %s", dimensions)

            calcStats = {}

            expSuccess = expStats.pop('SUCCESS')
            if expSuccess:
                calcStats['SUCCESS'] = True
            else:
                calcStats['FAILURES'] = list()
                success = False
                for detName, detStats in expStats.items():
                    detSuccess = detStats.pop('SUCCESS')
                    if not detSuccess:
                        for testName in expStats[detName]['FAILURES']:
                            calcStats['FAILURES'].append(detName + " " + testName)

            outputStats[expId] = calcStats

        runSuccess = True
        if len(self.config.runStatKeywords):
            outputStats['VERIFY'], runSuccess = self.verify(outputStats)

        outputStats['SUCCESS'] = success & runSuccess

        outputResults = mergeTable(inputResults, outputStats)
        outputMatrix = mergeTable(inputMatrix, outputStats)
        return pipeBase.Struct(
            outputStats=outputStats,
            outputResults=outputResults,
            outputMatrix=outputMatrix,
        )

    def verify(self, statisticsDictionary):
        """Verify if the measured statistics meet the verification criteria.

        Parameters
        ----------
        statisticsDictionary : `dict` [`str`, `dict`],
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
    camera = cT.PrerequisiteInput(
        name="camera",
        storageClass="Camera",
        doc="Input camera.",
        dimensions=["instrument", ],
        isCalibration=True,
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
        dimensions=["instrument", "exposure"],
        multiple=True,
    )
    inputMatrix = cT.Input(
        name="exposureMatrix",
        doc="Input matrix table to merge.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "exposure"],
        multiple=True,

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


class CpVerifyCalibMergeConfig(pipeBase.PipelineTaskConfig,
                               pipelineConnections=CpVerifyCalibMergeConnections):
    """Configuration paramters for exposure stats merging.
    """
    runStatKeywords = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Dictionary of statistics to run on the set of exposure values. The key should be the test "
        "name to record in the output, and the value should be the `lsst.afw.math` statistic name string.",
        default={},
    )
    hasMatrixCatalog = pexConfig.Field(
        dtype=bool,
        doc="Is there matrix catalog to merge?",
        default=False,
    )



class CpVerifyCalibMergeTask(pipeBase.PipelineTask):
    """Merge statistics from detectors together.
    """
    ConfigClass = CpVerifyCalibMergeConfig
    _DefaultName = 'cpVerifyCalibMerge'

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        dimensions = [dict(exp.dataId.required) for exp in inputRefs.inputStats]
        inputs['inputDims'] = dimensions

        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def run(self, inputStats, inputDims, inputResults=None, inputMatrix=None):
        """Merge statistics.

        Parameters
        ----------
        inputStats : `list` [`dict`]
            Measured statistics for a detector.
        inputDims : `list` [`dict`]
            List of dictionaries of input data dimensions/values.
            Each list entry should contain:

            ``"detector"``
                detector id value (`int`)

        Returns
        -------
        outputStats


        Notes
        -----
        The outputStats should have a yaml representation as follows.

        Detector detId:
          FAILURES:
          - Detector detId TEST_NAME
        SUCCESS: boolean
        """
        outputStats = {}
        success = True
        for detStats, dimensions in zip(inputStats, inputDims):
            detId = dimensions['detector']
            detName = f"Detector {detId}"
            calcStats = {}

            detSuccess = detStats.pop('SUCCESS')
            if detSuccess:
                calcStats['SUCCESS'] = True
            else:
                calcStats['FAILURES'] = list()
                success = False
                for testName in detStats['VERIFY']:
                    calcStats['FAILURES'].append(detName + " " + testName)

            outputStats[detName] = calcStats

        runSuccess = True
        if len(self.config.runStatKeywords):
            outputStats['VERIFY'], runSuccess = self.verify(outputStats)

        outputStats['SUCCESS'] = success & runSuccess

        return pipeBase.Struct(
            outputStats=outputStats,
        )

    def verify(self, statisticsDictionary):
        """Verify if the measured statistics meet the verification criteria.

        Parameters
        ----------
        statisticsDictionary : `dict` [`str`, `dict`],
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


def mergeTable(inputResults, newStats):
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
        return None

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
    if newStats:
        return vstack(vstack(inputResults), newStats)
    else:
        return vstack(inputResults)
