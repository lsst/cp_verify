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
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
import lsst.pex.config as pexConfig


__all__ = ['CpVerifyExpMergeConfig', 'CpVerifyExpMergeTask',
           'CpVerifyRunMergeConfig', 'CpVerifyRunMergeTask']


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


class CpVerifyExpMergeConfig(pipeBase.PipelineTaskConfig,
                             pipelineConnections=CpVerifyExpMergeConnections):
    """Configuration paramters for exposure stats merging.
    """
    exposureStatKeywords = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Dictionary of statistics to run on the set of detector values. The key should be the test "
        "name to record in the output, and the value should be the `lsst.afw.math` statistic name string."
        default={},
    )


class CpVerifyExpMergeTask(pipeBase.PipelineTask, pipeBase.CmdLineTask):
    """Merge statistics from detectors together.
    """
    ConfigClass = CpVerifyExpMergeConfig
    _DefaultName = 'cpVerifyExpMerge'

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        dimensions = [exp.dataId.byName() for exp in inputRefs.inputStats]
        inputs['inputDims'] = dimensions

        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def run(self, inputStats, camera, inputDims):
        """Merge statistics.

        Parameters
        ----------
        inputStats : `list` [`dict`]
            Measured statistics for a detector.
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
        outputStats : `dict`
            Merged full exposure statistics.

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

        for detStats, dimensions in zip(inputStats, inputDims):
            detId = dimensions['detector']
            detName = camera[detId].getName()
            calcStats = {}

            if detStats['SUCCESS'] is True:
                calcStats['SUCCESS'] = True
            else:
                calcStats['SUCCESS'] = False
                calcStats['FAILURES'] = list()
                success = False
                # See if the detector failed
                if 'DET' in detStats['VERIFY']:
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
                        for testName, testResult in ampStats:
                            if testResult is False:
                                calcStats['FAILURES'].append(ampName + " " + testName)

            outputStats[detName] = calcStats

        exposureSuccess = True
        if len(self.config.exposureStatKeywords):
            outputStats['VERIFY'], exposureSuccess = self.verify(outputStats)

        outputStats['SUCCESS'] = success & exposureSuccess

        return pipeBase.Struct(
            outputStats=outputStats,
        )

    def verify(self, statisticsDictionary):
        """Verify if the measured statistics meet the verification criteria.

        Parameters
        ----------
        statisticsDictionary : `dict` [`str`, `dict` [`str`, scalar]],
            Dictionary of measured statistics.  The inner dictionary
            should have keys that are statistic names (`str`) with
            values that are some sort of scalar (`int` or `float` are
            the mostly likely types).

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

    outputStats = cT.Output(
        name="runStats",
        doc="Output statistics.",
        storageClass="StructuredDataDict",
        dimensions=["instrument"],
    )


class CpVerifyRunMergeConfig(pipeBase.PipelineTaskConfig,
                             pipelineConnections=CpVerifyRunMergeConnections):
    """Configuration paramters for exposure stats merging.
    """
    runStatKeywords = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Dictionary of statistics to run on the set of exposure values. The key should be the test "
        "name to record in the output, and the value should be the `lsst.afw.math` statistic name string."
        default={},
    )


class CpVerifyRunMergeTask(pipeBase.PipelineTask, pipeBase.CmdLineTask):
    """Merge statistics from detectors together.
    """
    ConfigClass = CpVerifyRunMergeConfig
    _DefaultName = 'cpVerifyRunMerge'

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        dimensions = [exp.dataId.byName() for exp in inputRefs.inputStats]
        inputs['inputDims'] = dimensions

        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def run(self, inputStats, inputDims):
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

        Returns
        -------
        outputStats : `dict`
            Merged full exposure statistics.

        Notes
        -----
        The outputStats should have a yaml representation as follows.

        VERIFY:
          ExposureId1:
            VERIFY_MEAN: boolean
            VERIFY_SIGMA: boolean
          ExposureId2:
            [...]
          MEAN_UNIMODAL: boolean
          SIGMA_UNIMODAL: boolean
        """
        outputStats = {}
        success = True
        for expStats, dimensions in zip(inputStats, inputDims):
            expId = dimensions['exposure']
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
