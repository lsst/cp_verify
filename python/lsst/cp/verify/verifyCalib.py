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
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT

from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration

__all__ = ['CpVerifyCalibConfig', 'CpVerifyCalibTask']


class CpVerifyCalibConnections(pipeBase.PipelineTaskConnections,
                               dimensions={"instrument", "detector"},
                               defaultTemplates={}):
    inputCalib = cT.Input(
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
        lookupFunction=lookupStaticCalibration,
        isCalibration=True,
    )

    outputStats = cT.Output(
        name="calibStats",
        doc="Output statistics from cp_verify.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "detector"],
    )


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


class CpVerifyCalibTask(pipeBase.PipelineTask):
    """Main statistic measurement and validation class.

    This operates on a generic calibration, and is designed to be
    subclassed so specific calibrations can apply their own validation
    methods.
    """

    ConfigClass = CpVerifyCalibConfig
    _DefaultName = 'cpVerifyCalib'

    def run(self, inputCalib, camera=None):
        """Calculate quality statistics and verify they meet the requirements
        for a calibration.

        Parameters
        ----------
        inputCalib : `lsst.ip.isr.IsrCalib`
            The calibration to be measured.
        camera : `lsst.afw.cameraGeom.Camera`, optional
            Input camera.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:
            - ``outputStats`` : `dict`
                The output measured statistics.

        Notes
        -----
        The outputStats should have a yaml representation of the form
        (with STAT and TEST being the appropriate statistic and test
        names)

        DET:
          STAT: value
          STAT2: value
        AMP:
          STAT: value
          STAT2: value
        VERIFY:
          TEST: boolean
        SUCCESS: boolean

        """
        outputStats = {}
        outputStats['AMP'] = self.amplifierStatistics(inputCalib, camera=camera)
        outputStats['DET'] = self.detectorStatistics(inputCalib, camera=camera)
        outputStats['VERIFY'], outputStats['SUCCESS'] = self.verify(inputCalib, outputStats, camera=camera)

        return pipeBase.Struct(
            outputStats=outputStats,
        )

    # Methods that need to be implemented by the calibration-level subclasses.
    def detectorStatistics(self, inputCalib, camera=None):
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

        Raises
        ------
        NotImplementedError :
            This method must be implemented by the calibration-type
            subclass.
        """
        raise NotImplementedError("Subclasses must implement detector statistics method.")

    def amplifierStatistics(self, inputCalib, camera=None):
        """Calculate amplifier level statistics from the calibration.

        Parameters
        ----------
        inputCalib : `lsst.ip.isr.IsrCalib`
            The calibration to verify.
        camera : `lsst.afw.cameraGeom.Camera`, optional
            Input camera.

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

    def verify(self, inputCalib, statisticsDict, camera=None):
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
