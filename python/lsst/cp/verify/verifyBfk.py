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
import lsst.pipe.base.connectionTypes as cT

from lsstDebug import getDebugFrame
from lsst.cp.pipe.utils import (funcPolynomial, irlsFit)
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

    def setDefaults(self):
        super().setDefaults()
        self.catalogStatKeywords = {'BRIGHT_SLOPE': None}  # noqa F841


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
        magnitude = -2.5 * np.log10(catalog.getPsfInstFlux())
        # This is a simple version
        # size = np.sqrt(0.5 * (catalog.getIxx() + catalog.getIyy()))
        # This is the pipe_analysis version
        srcSize = np.sqrt(0.5*(catalog['base_SdssShape_xx'] + catalog['base_SdssShape_yy']))
        psfSize = np.sqrt(0.5*(catalog['base_SdssShape_psf_xx'] + catalog['base_SdssShape_psf_yy']))
        size = 100*(srcSize - psfSize)/(0.5*(srcSize + psfSize))

        if 'BRIGHT_SLOPE' in self.config.catalogStatKeywords:
            linearFit, linearFitErr, chiSq, weights = irlsFit([np.median(size), 0.0],
                                                              magnitude, size, funcPolynomial)

            outputStatistics['BRIGHT_SLOPE'] = float(linearFit[1])
            self.debugFit('brightSlope', magnitude, size, linearFit)

        # This should have a debug view
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
            #            import pdb; pdb.set_trace()
            figure, axes = plt.subplots(1)

            axes.scatter(xVector, yVector, c='blue', marker='+')
            modelX = np.arange(np.min(xVector) - 1, np.max(xVector) + 1, 0.1)
            axes.plot(modelX, yModel[0] + yModel[1] * modelX, 'r-')
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

        catalogVerify['BRIGHT_SLOPE'] = True
        if np.abs(catStats['BRIGHT_SLOPE']) > 0.05:
            catalogVerify['BRIGHT_SLOPE'] = False

        return {'AMP': verifyStats, 'CATALOG': catalogVerify}, bool(success)
