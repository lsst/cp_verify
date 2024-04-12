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

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT

__all__ = ['CpVerifyCalibConfig', 'CpVerifyCalibTask']


class CpVerifyCalibConnections(pipeBase.PipelineTaskConnections,
                               dimensions={"instrument", "detector"},
                               defaultTemplates={}):

    exposure = cT.Input(
        name="raw",
        doc="Exposure to retreve calibration",
        storageClass='Exposure',
        dimensions=("instrument", "detector", "exposure"),
        multiple=True,
        deferLoad=True,
    )

    inputCalib = cT.PrerequisiteInput(
        name="calib",
        doc="Input calib to calculate statistics for.",
        storageClass="IsrCalib",
        dimensions=["instrument", "detector"],
        isCalibration=True
    )

    camera = cT.PrerequisiteInput(
        name="camera",
        doc="Input camera to use for gain lookup.",
        storageClass="Camera",
        dimensions=("instrument",),
        isCalibration=True,
    )

    outputStats = cT.Output(
        name="calibStats",
        doc="Output statistics from cp_verify.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "detector"],
    )
    outputResults = cT.Output(
        name="detectorResults",
        doc="Output results from cp_verify.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "detector"],
    )
    outputMatrix = cT.Output(
        name="detectorMatrix",
        doc="Output matrix results from cp_verify.",
        storageClass="ArrowAstropy",
        dimensions=["instrument", "detector"],
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not config.hasMatrixCatalog:
            self.outputs.discard("outputMatrix")


class CpVerifyCalibConfig(pipeBase.PipelineTaskConfig,
                          pipelineConnections=CpVerifyCalibConnections):
    """Configuration parameters for CpVerifyCalibTask.
    """
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
    calibStatKeywords = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        doc="Calib statistics to run.",
        default={},
    )

    stageName = pexConfig.Field(
        dtype=str,
        doc="Stage name to use in columns.",
        default="NOCALIB",
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


class CpVerifyCalibTask(pipeBase.PipelineTask):
    """Main statistic measurement and validation class.

    This operates on a generic calibration, and is designed to be
    subclassed so specific calibrations can apply their own validation
    methods.
    """

    ConfigClass = CpVerifyCalibConfig
    _DefaultName = 'cpVerifyCalib'

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        inputs["dimensions"] = dict(inputRefs.inputCalib.dataId.required)

        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def run(self,
            inputCalib,
            camera=None,
            exposure=None,
            dimensions=None,
            ):
        """Calculate quality statistics and verify they meet the requirements
        for a calibration.

        Parameters
        ----------
        inputCalib : `lsst.ip.isr.IsrCalib`
            The calibration to be measured.
        camera : `lsst.afw.cameraGeom.Camera`, optional
            Input camera.
        exposure : `lsst.afw.image.Exposure`, optional
            Dummy exposure to identify a particular calibration
            dataset.
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
        outputStats['AMP'] = self.amplifierStatistics(inputCalib, camera=camera)
        outputStats['DET'] = self.detectorStatistics(inputCalib, camera=camera)
        outputStats['VERIFY'], outputStats['SUCCESS'] = self.verify(inputCalib, outputStats, camera=camera)

        outputResults, outputMatrix = self.repackStats(outputStats, dimensions)
        if outputResults is not None:
            outputResults = Table(outputResults)
        if outputMatrix is not None:
            outputMatrix = Table(outputMatrix)

        return pipeBase.Struct(
            outputStats=outputStats,
            outputResults=outputResults,
            outputMatrix=outputMatrix,
        )

    # Methods that need to be implemented by the calibration-level subclasses.
    def detectorStatistics(self, inputCalib, camera=None, exposure=None):
        """Calculate detector level statistics from the calibration.

        Parameters
        ----------
        inputCalib : `lsst.ip.isr.IsrCalib`
            The calibration to verify.

        Returns
        -------
        outputStatistics : `dict` [`str`, scalar]
            A dictionary of the statistics measured and their values.
        camera : `lsst.afw.cameraGeom.Camera`, optional
            Input camera.
        exposure : `lsst.afw.image.Exposure`, optional
            Dummy exposure to identify a particular calibration
            dataset.

        Raises
        ------
        NotImplementedError :
            This method must be implemented by the calibration-type
            subclass.
        """
        raise NotImplementedError("Subclasses must implement detector statistics method.")

    def amplifierStatistics(self, inputCalib, camera=None, exposure=None):
        """Calculate amplifier level statistics from the calibration.

        Parameters
        ----------
        inputCalib : `lsst.ip.isr.IsrCalib`
            The calibration to verify.
        camera : `lsst.afw.cameraGeom.Camera`, optional
            Input camera.
        exposure : `lsst.afw.image.Exposure`, optional
            Dummy exposure to identify a particular calibration
            dataset.

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
        raise NotImplementedError("Subclasses must implement amplifier statistics method.")

    def verify(self, inputCalib, statisticsDict, camera=None, exposure=None):
        """Verify that the measured calibration meet the verification criteria.

        Parameters
        ----------
        inputCalib : `lsst.ip.isr.IsrCalib`
            The calibration to verify.
        statisticsDictionary : `dict` [`str`, `dict` [`str`, scalar]],
            Dictionary of measured statistics.  The inner dictionary
            should have keys that are statistic names (`str`) with
            values that are some sort of scalar (`int` or `float` are
            the mostly likely types).
        camera : `lsst.afw.cameraGeom.Camera`, optional
            Input camera.
        exposure : `lsst.afw.image.Exposure`, optional
            Dummy exposure to identify a particular calibration
            dataset.

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

        This method should be redefined in subclasses.

        Parameters
        ----------
        statisticsDictionary : `dict` [`str`, `dict` [`str`, scalar]],
            Dictionary of measured statistics.  The inner dictionary
            should have keys that are statistic names (`str`) with
            values that are some sort of scalar (`int` or `float` are
            the mostly likely types).

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
        for ampName, stats in statisticsDict["VERIFY"]["AMP"].items():
            for key, value in stats.items():
                rows[ampName][f"{self.config.stageName}_VERIFY_{key}"] = value

        # pack final list
        for ampName, stats in rows.items():
            rowList.append(stats)

        return rowList, matrixRowList
