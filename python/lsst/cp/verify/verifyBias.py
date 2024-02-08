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

import lsst.afw.math as afwMath

from astropy.table import Table

from lsst.geom import Point2I, Extent2I, Box2I
from lsst.pex.config import Field
from .verifyStats import CpVerifyStatsConfig, CpVerifyStatsTask, CpVerifyStatsConnections

__all__ = ['CpVerifyBiasConfig', 'CpVerifyBiasTask']


class CpVerifyBiasConfig(CpVerifyStatsConfig,
                         pipelineConnections=CpVerifyStatsConnections):
    """Inherits from base CpVerifyStatsConfig.
    """

    ampCornerBoxSize = Field(
        dtype=int,
        doc="Size of box to use for measure corner signal.",
        default=200,
    )

    def setDefaults(self):
        super().setDefaults()
        self.imageStatKeywords = {'MEAN': 'MEAN',  # noqa F841
                                  'NOISE': 'STDEVCLIP', }
        self.crImageStatKeywords = {'CR_NOISE': 'STDEV', }  # noqa F841
        self.metadataStatKeywords = {'RESIDUAL STDEV': 'AMP', }  # noqa F841


class CpVerifyBiasTask(CpVerifyStatsTask):
    """Bias verification sub-class, implementing the verify method.
    """
    ConfigClass = CpVerifyBiasConfig
    _DefaultName = 'cpVerifyBias'

    stageName = "bias"

    def imageStatistics(self, exposure, uncorrectedExposure, statControl):
        # Docstring inherited
        outputStatistics = super().imageStatistics(exposure, uncorrectedExposure, statControl)

        boxSize = self.config.ampCornerBoxSize
        statisticToRun = afwMath.stringToStatisticsProperty("MEAN")

        for ampIdx, amp in enumerate(exposure.getDetector()):
            ampName = amp.getName()

            bbox = amp.getBBox()
            xmin = bbox.getMaxX() - boxSize if amp.getRawFlipX() else bbox.getMinX()
            ymin = bbox.getMaxY() - boxSize if amp.getRawFlipY() else bbox.getMinY()
            llc = Point2I(xmin, ymin)
            extent = Extent2I(boxSize, boxSize)
            cornerBox = Box2I(llc, extent)
            cornerExp = exposure[cornerBox]

            stats = afwMath.makeStatistics(
                cornerExp.getMaskedImage(), statisticToRun, statControl
            )
            outputStatistics[ampName]['AMP_CORNER'] = stats.getValue()

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
        detector = exposure.getDetector()
        ampStats = statisticsDict['AMP']
        metadataStats = statisticsDict['METADATA']

        verifyStats = {}
        success = True
        for ampName, stats in ampStats.items():
            verify = {}

            # DMTN-101 Test 4.2: Mean is 0.0 within noise.
            verify['MEAN'] = bool(np.abs(stats['MEAN']) < stats['NOISE'])

            # DMTN-101 Test 4.3: Clipped mean matches readNoise.  This
            # test should use the nominal detector read noise.  The
            # f"RESIDUAL STDEV {ampName}" metadata entry contains the
            # measured dispersion in the overscan-corrected overscan
            # region, which should provide an estimate of the read
            # noise.  However, directly using this value will cause
            # some fraction of verification runs to fail if the
            # scatter in read noise values is comparable to the test
            # threshold, as the overscan residual measured may be
            # sampling from the low end tail of the distribution.
            # This measurement is also likely to be smaller than that
            # measured on the bulk of the image as the overscan
            # correction should be an optimal fit to the overscan
            # region, but not necessarily for the image region.
            readNoise = detector[ampName].getReadNoise()
            verify['NOISE'] = bool((stats['NOISE'] - readNoise)/readNoise <= 0.05)

            # DMTN-101 Test 4.4: CR rejection matches clipped mean.
            verify['CR_NOISE'] = bool(np.abs(stats['NOISE'] - stats['CR_NOISE'])/stats['CR_NOISE'] <= 0.05)

            # Confirm this hasn't triggered a raise condition.
            if 'FORCE_FAILURE' in stats:
                verify['PROCESSING'] = False

            verify['SUCCESS'] = bool(np.all(list(verify.values())))
            if verify['SUCCESS'] is False:
                success = False

            # After determining the verification status for this
            # exposure, we can also check to see how well the read
            # noise measured from the overscan residual matches the
            # nominal value used above in Test 4.3.  If these disagree
            # consistently and significantly, then the assumptions
            # used in that test may be incorrect, and the nominal read
            # noise may need recalculation.  Only perform this check
            # if the metadataStats contain the required entry.
            if 'RESIDUAL STDEV' in metadataStats and ampName in metadataStats['RESIDUAL STDEV']:
                verify['READ_NOISE_CONSISTENT'] = True
                overscanReadNoise = metadataStats['RESIDUAL STDEV'][ampName]
                if overscanReadNoise:
                    if ((overscanReadNoise - readNoise)/readNoise > 0.05):
                        verify['READ_NOISE_CONSISTENT'] = False

            verifyStats[ampName] = verify

        return {'AMP': verifyStats}, bool(success)

    def repackStats(self, detStats, detDims):
        """Repack hierarchical results into flat table.

        Parameters
        ----------
        detStats : `dict` [`str`, `dict`]
            A nested set of dictionaries containing relevant
            statistics.
        detDims : `dict` [`str`, `str`]
            The dimensions of this set of statistics.

        Returns
        -------
        outputResults : `astropy.Table`, optional
            The repacked flat table.
        outputMatrix : `astropy.Table`, optional
            The repackaed matrix data, in a flat table.

        Raises
        ------
        NotImplementedError :
            This method must be implemented by the calibration-type
            subclass.
        """
        rowList = []

        row = {}
        instrument = detDims["instrument"]
        exposure = detDims["exposure"]
        detector = detDims["detector"]
        mjd = detStats["ISR"]["MJD"] if "ISR" in detStats and "MJD" in detStats["ISR"] else 0.0

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
        if "METADATA" in detStats:
            if "RESIDUAL STDEV" in detStats["METADATA"]:
                for ampName, value in detStats["METADATA"]["RESIDUAL STDEV"].items():
                    row[ampName]["biasReadNoise"] = value

        # Get verify stats
        for ampName, stats in detStats["VERIFY"]["AMP"].items():
            row[ampName]["biasVerifyMean"] = stats["MEAN"]
            row[ampName]["biasVerifyNoise"] = stats["NOISE"]
            row[ampName]["biasVerifyCrNoise"] = stats["CR_NOISE"]
            row[ampName]["biasVerifyReadNoiseConsistent"] = stats.get("READ_NOISE_CONSISTENT", False)

        # Get isr stats
        if "ISR" in detStats:
            for ampName, stats in detStats["ISR"]["CALIBDIST"].items():
                for level in self.config.expectedDistributionLevels:
                    key = f"LSST CALIB {self.stageName.upper()} {ampName} DISTRIBUTION {level}-PCT"
                    row[ampName][f"biasDistribution_{level}"] = stats[key]

            projStats = detStats["ISR"]["PROJECTION"]
            for ampName in projStats["AMP_HPROJECTION"].keys():
                row[ampName]["biasSerialProfile"] = np.array(projStats["AMP_HPROJECTION"][ampName])
            for ampName in projStats["AMP_VPROJECTION"].keys():
                row[ampName]["biasParallelProfile"] = np.array(projStats["AMP_VPROJECTION"][ampName])

        # Create output table:
        for ampName, stats in row.items():
            rowList.append(stats)
        return Table(rowList), None

        # # We need all rows of biasParallelProfile and biasParallelProfile
        # # to be the same length for serialization. Therefore, we pad
        # # to the longest length.

        # maxSerialLen = 0
        # maxParallelLen = 0

        # for row in rowList:
        #     if len(row["biasSerialProfile"]) > maxSerialLen:
        #         maxSerialLen = len(row["biasSerialProfile"])
        #     if len(row["biasParallelProfile"]) > maxParallelLen:
        #         maxParallelLen = len(row["biasParallelProfile"])

        # for row in rowList:
        #     if len(row["biasSerialProfile"]) < maxSerialLen:
        #         row["biasSerialProfile"] = np.pad(
        #             row["biasSerialProfile"],
        #             (0, maxSerialLen - len(row["biasSerialProfile"])),
        #             constant_values=np.nan,
        #         )
        #     if len(row["biasParallelProfile"]) < maxParallelLen:
        #         row["biasParallelProfile"] = np.pad(
        #             row["biasParallelProfile"],
        #             (0, maxParallelLen - len(row["biasParallelProfile"])),
        #             constant_values=np.nan,
        #         )
