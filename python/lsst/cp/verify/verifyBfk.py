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
import lsst.afw.table as afwTable
import lsst.pex.config as pexConfig
import lsst.pipe.base.connectionTypes as cT

from scipy.optimize import least_squares
from lsstDebug import getDebugFrame
from .verifyStats import CpVerifyStatsConfig, CpVerifyStatsTask, CpVerifyStatsConnections


__all__ = ['CpVerifyBfkConfig', 'CpVerifyBfkTask']


class CpVerifyBfkConnections(CpVerifyStatsConnections,
                             dimensions={'instrument', 'visit', 'detector'}):
    inputExp = cT.Input(
        name="icExp",
        doc="Input exposure to calculate statistics for.",
        storageClass="ExposureF",
        dimensions=["instrument", "visit", "detector"],
    )
    inputCatalog = cT.Input(
        name="icSrc",
        doc="Input catalog to calculate statistics from.",
        storageClass="SourceCatalog",
        dimensions=["instrument", "visit", "detector"],
    )
    uncorrectedCatalog = cT.Input(
        name="uncorrectedSrc",
        doc="Input catalog without correction applied.",
        storageClass="SourceCatalog",
        dimensions=["instrument", "visit", "detector"],
    )
    camera = cT.PrerequisiteInput(
        name="camera",
        storageClass="Camera",
        doc="Input camera.",
        dimensions=["instrument", ],
        isCalibration=True,
    )

    outputStats = cT.Output(
        name="detectorStats",
        doc="Output statistics from cp_verify.",
        storageClass="StructuredDataDict",
        dimensions=["instrument", "visit", "detector"],
    )


class CpVerifyBfkConfig(CpVerifyStatsConfig,
                        pipelineConnections=CpVerifyBfkConnections):
    """Inherits from base CpVerifyStatsConfig.
    """

    matchRadiusPix = pexConfig.Field(
        dtype=float,
        default=3,
        doc=("Match radius for matching icSourceCat objects to sourceCat objects (pixels)"),
    )

    def setDefaults(self):
        super().setDefaults()
        self.stageName = 'BFK'
        self.catalogStatKeywords = {'BRIGHT_SLOPE': None,
                                    'NUM_MATCHES': None,
                                   }  # noqa F841


def exponentialModel(x, A, alpha, x0):
    """An exponential model.
    """
    return A * np.exp(alpha * (x - x0))


def modelResidual(p, x, y):
    """Model residual for fit below.
    """
    return y - exponentialModel(x, *p)


class CpVerifyBfkTask(CpVerifyStatsTask):
    """Bfk verification sub-class, implementing the verify method.
    """

    ConfigClass = CpVerifyBfkConfig
    _DefaultName = 'cpVerifyBfk'

    def catalogStatistics(self, exposure, catalog, uncorrectedCatalog, statControl):
        """Measure the catalog statistics.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            The exposure to measure.
        catalog : `lsst.afw.table.Table`
            The catalog to measure.
        uncorrectedCatalog : `lsst.afw.table.Table`
            The uncorrected catalog to measure.
        statControl : `lsst.afw.math.StatisticsControl`
            Statistics control object with parameters defined by
            the config.

        Returns
        -------
        outputStatistics : `dict` [`str`, `dict` [`str`, scalar]]
            A dictionary indexed by the amplifier name, containing
            dictionaries of the statistics measured and their values.
        """
        outputStatistics = {}

        matches = afwTable.matchXy(catalog, uncorrectedCatalog, self.config.matchRadiusPix)
        outputStatistics['NUM_MATCHES'] = len(matches)

        magnitude = []
        sizeDiff = []
        if not matches:
            outputStatistics['MAGNITUDES'] = magnitude
            outputStatistics['SIZE_DIFF'] = sizeDiff
            return outputStatistics

        xxKey = 'ext_shapeHSM_HsmSourceMoments_xx'
        yyKey = 'ext_shapeHSM_HsmSourceMoments_yy'
        for source, uncorrectedSource, d in matches:
            # This uses the simple difference in source moments.
            sourceMagnitude = -2.5 * np.log10(source.getPsfInstFlux())
            sourceSize = source[xxKey] + source[yyKey]
            uncorrectedSize = uncorrectedSource[xxKey] + uncorrectedSource[yyKey]

            magnitude.append(float(sourceMagnitude))
            sizeDiff.append(float(uncorrectedSize - sourceSize))

        mask = np.isfinite(magnitude) * np.isfinite(sizeDiff)
        if 'BRIGHT_SLOPE' in self.config.catalogStatKeywords:
            exponentialFit = least_squares(modelResidual, [0.0, -0.01, 0.0],
                                           args=(np.array(magnitude)[mask], np.array(sizeDiff)[mask]),
                                           loss='cauchy')

            outputStatistics['BRIGHT_SLOPE'] = float(exponentialFit.x[1])
            outputStatistics['FIT_SUCCESS'] = exponentialFit.success
            self.debugFit('brightSlope', magnitude, sizeDiff, exponentialFit.x)

        outputStatistics['MAGNITUDES'] = magnitude
        outputStatistics['SIZE_DIFF'] = sizeDiff

        return outputStatistics

    def debugFit(self, stepname, xVector, yVector, yModel):
        """Debug method for linearity fitting.
        Parameters
        ----------
        stepname : `str`
            A label to use to check if we care to debug at a given
            line of code.
        xVector : `numpy.array`, (N,)
            The values to use as the independent variable in the
            linearity fit.
        yVector : `numpy.array`, (N,)
            The values to use as the dependent variable in the
            linearity fit.
        yModel : `numpy.array`, (N,)
            The values to use as the linearized result.
        """
        frame = getDebugFrame(self._display, stepname)
        if frame:
            import matplotlib.pyplot as plt
            figure, axes = plt.subplots(1)

            axes.scatter(xVector, yVector, c='blue', marker='+')
            modelX = np.arange(np.min(xVector) - 1, np.max(xVector) + 1, 0.1)
            axes.plot(modelX, exponentialModel(modelX, *yModel), 'r-')
            plt.xlabel("Instrumental PSF magnitude")
            plt.ylabel("Source size trace")
            plt.title(f"BFK slope: {yModel[0]:.3f} + {yModel[1]:.3f} m")
            figure.show()
            prompt = "Press Enter or c to continue [chp]..."

            while True:
                ans = input(prompt).lower()
                if ans in ("", " ", "c",):
                    break
                elif ans in ("p", ):
                    import pdb
                    pdb.set_trace()
                elif ans in ("h", ):
                    print("[h]elp [c]ontinue [p]db e[x]itDebug")
            plt.close()

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
        verifyStats = {}
        success = True
        catalogVerify = {}
        catStats = statisticsDict['CATALOG']

        catalogVerify['BRIGHT_SLOPE'] = False
        catalogVerify['NUM_MATCHES'] = False
        # These values need justification.
        if catStats['FIT_SUCCESS'] and catStats['BRIGHT_SLOPE'] < -0.5:
            catalogVerify['BRIGHT_SLOPE'] = True
        if catStats['NUM_MATCHES'] > 10:
            catalogVerify['NUM_MATCHES'] = True

        if catalogVerify['NUM_MATCHES'] and not catalogVerify['BRIGHT_SLOPE']:
            success = False
        return {'AMP': verifyStats, 'CATALOG': catalogVerify}, bool(success)

    def repackStats(self, statisticsDict, dimensions):
        # docstring inherited
        rows = {}
        rowList = []
        matrixRowList = None

        rowBase = {
            "instrument": dimensions["instrument"],
            "detector": dimensions["detector"],
        }

        # CAT results
        for ampName, stats in statisticsDict["CAT"].items():
            rows[ampName] = {}
            rows[ampName].update(rowBase)
            rows[ampName]["amplifier"] = ampName
            for key, value in stats.items():
                rows["ampName"][f"{self.config.stageName}_{key}"] = value
                rows["ampName"][f"{self.config.stageName}_VERIFY_{key}"] = value

        # pack final list
        for ampName, stats in rows.items():
            rowList.append(stats)

        return rowList, matrixRowList
