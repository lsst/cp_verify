# This file is part of cp_verify.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import unittest

import lsst.utils.tests
import lsst.ip.isr.isrMock as isrMock
import lsst.cp.verify as cpVerify
import lsst.ip.isr.isrFunctions as isrFunctions


def updateMockExp(exposure, addCR=True):
    """Update an exposure with a mask and variance plane.

    Parameters
    ----------
    exposure : `lsst.afw.image.Exposure`
        Exposure to be modified in place.
    addCR : `bool`
        Whether a known cosmic ray should be added to ``exposure``.
    """
    if addCR:
        # Add a cosmic ray
        image = exposure.getImage()
        image.getArray()[50, 50] = 10000.0

    # Set the mask and variance planes:
    mask = exposure.getMask()
    mask.getArray()[:, 10] = 1
    isrFunctions.updateVariance(exposure.getMaskedImage(), 1.0, 5.0)


class ToySubClass(cpVerify.CpVerifyStatsTask):
    """The CpVerifyStatsTask requires an implentation of verify.
    """

    def verify(self, inputExp, outputStats):
        # Docstring inherited from CpVerifyStatsTask.verify()
        verifiedStats = {'A REAL TEST': True, 'A BAD TEST': False}
        successValue = True

        return verifiedStats, successValue


class VerifyStatsTestCase(lsst.utils.tests.TestCase):
    """Unit test for stats code.
    """

    def setUp(self):
        """Generate a mock exposure/camera to test."""
        self.inputExp = isrMock.CalibratedRawMock().run()
        self.camera = isrMock.IsrMock().getCamera()

        updateMockExp(self.inputExp)

    def test_failures(self):
        """Test that all the NotImplementedError methods fail correctly."""
        results = None
        with self.assertRaises(NotImplementedError):
            # We have not implemented a verify method
            config = cpVerify.CpVerifyStatsConfig()
            config.numSigmaClip = 3.0
            task = cpVerify.CpVerifyStatsTask(config=config)
            results = task.run(self.inputExp, self.camera)

            # Or the catalog stats
            config.catalogStatKeywords = {'CAT_MEAN', 'MEDIAN'}
            task = cpVerify.CpVerifyStatsTask(config=config)
            results = task.run(self.inputExp, self.camera)

            # Or the detector stats
            config.catalogStatKeywords = {}
            config.detectorStatKeywords = {'DET_SIGMA', 'STDEV'}
            task = cpVerify.CpVerifyStatsTask(config=config)
            results = task.run(self.inputExp, self.camera)
        self.assertIsNone(results)

    def test_generic(self):
        """Test a subset of the output values to identify that the
        image stat methods haven't changed.
        """
        config = cpVerify.CpVerifyStatsConfig()
        config.imageStatKeywords = {'MEAN': 'MEAN', 'MEDIAN': 'MEDIAN', 'CLIPPED': 'MEANCLIP',
                                    'SIGMA': 'STDEV'}
        config.unmaskedImageStatKeywords = {'un_MEAN': 'MEAN', 'un_MEDIAN': 'MEDIAN',
                                            'un_CLIPPED': 'MEANCLIP',
                                            'un_SIGMA': 'STDEV'}
        config.crImageStatKeywords = {'cr_MEAN': 'MEAN', 'cr_MEDIAN': 'MEDIAN', 'cr_CLIPPED': 'MEANCLIP',
                                      'cr_SIGMA': 'STDEV'}
        config.normImageStatKeywords = {'norm_MEAN': 'MEAN', 'norm_MEDIAN': 'MEDIAN',
                                        'norm_CLIPPED': 'MEANCLIP',
                                        'norm_SIGMA': 'STDEV'}
        config.numSigmaClip = 3.0
        task = ToySubClass(config=config)

        results = task.run(self.inputExp, self.camera)
        resultStats = results.outputStats

        self.assertAlmostEqual(resultStats['AMP']['C:0,0']['MEAN'], 1506.06976, 4)
        self.assertAlmostEqual(resultStats['AMP']['C:0,0']['un_MEAN'], 1501.0299, 4)
        self.assertAlmostEqual(resultStats['AMP']['C:0,0']['norm_MEAN'], 301.213957, 4)
        self.assertAlmostEqual(resultStats['AMP']['C:0,0']['cr_MEAN'], 1504.2776, 4)

        self.assertTrue(resultStats['VERIFY']['A REAL TEST'])
        self.assertFalse(resultStats['VERIFY']['A BAD TEST'])

        self.assertTrue(resultStats['SUCCESS'])


class VerifyBiasTestCase(lsst.utils.tests.TestCase):
    """Unit test for stats code - bias cases."""

    def setUp(self):
        """Generate a mock exposure/camera to test."""
        config = isrMock.IsrMockConfig()
        config.isTrimmed = True
        config.rngSeed = 12345
        biasExposure = isrMock.BiasMock(config=config).run()

        config.rngSeed = 54321
        fakeBias = isrMock.BiasMock(config=config).run()

        self.inputExp = biasExposure.clone()
        mi = self.inputExp.getMaskedImage()
        mi.scaledMinus(1.0, fakeBias.getMaskedImage())
        updateMockExp(self.inputExp)

        self.camera = isrMock.IsrMock().getCamera()

    def test_bias(self):
        """Test a subset of the output values to identify that the
        image stat methods haven't changed.
        """
        config = cpVerify.CpVerifyBiasConfig()
        config.numSigmaClip = 3.0
        task = cpVerify.CpVerifyBiasTask(config=config)
        results = task.run(self.inputExp, self.camera)
        biasStats = results.outputStats

        self.assertAlmostEqual(biasStats['AMP']['C:0,0']['MEAN'], 2.08672, 4)
        self.assertAlmostEqual(biasStats['AMP']['C:0,0']['NOISE'], 13.99547, 4)
        self.assertAlmostEqual(biasStats['AMP']['C:0,0']['CR_NOISE'], 14.11526, 4)

        self.assertIn(biasStats['SUCCESS'], [True, False])


class VerifyDarkTestCase(lsst.utils.tests.TestCase):
    """Unit test for stats code - dark cases.
    """

    def setUp(self):
        """Generate a mock exposure/camera to test."""
        config = isrMock.IsrMockConfig()
        config.isTrimmed = True
        config.rngSeed = 12345
        darkExposure = isrMock.DarkMock(config=config).run()

        config.rngSeed = 54321
        fakeDark = isrMock.DarkMock(config=config).run()

        self.inputExp = darkExposure.clone()
        mi = self.inputExp.getMaskedImage()
        mi.scaledMinus(1.0, fakeDark.getMaskedImage())
        updateMockExp(self.inputExp)

        # Use this to test the metadataStats code, as this is the case
        # it's designed to fix.
        metadataContents = {}
        metadataContents["RESIDUAL STDEV C:0,0"] = 12.0
        metadataContents["RESIDUAL STDEV"] = 24.0
        self.metadata = {}
        self.metadata["subGroup"] = metadataContents

        self.camera = isrMock.IsrMock().getCamera()

    def test_dark(self):
        """Test a subset of the output values to identify that the
        image stat methods haven't changed.
        """
        config = cpVerify.CpVerifyDarkConfig()
        config.numSigmaClip = 3.0
        task = cpVerify.CpVerifyDarkTask(config=config)
        results = task.run(self.inputExp, self.camera, taskMetadata=self.metadata)
        darkStats = results.outputStats

        self.assertAlmostEqual(darkStats['AMP']['C:0,0']['MEAN'], 2.0043, 4)
        self.assertAlmostEqual(darkStats['AMP']['C:0,0']['NOISE'], 3.12948, 4)
        self.assertAlmostEqual(darkStats['AMP']['C:0,0']['CR_NOISE'], 3.15946, 4)

        self.assertIn(darkStats['SUCCESS'], [True, False])


class VerifyDefectsTestCase(lsst.utils.tests.TestCase):
    """Unit test for stats code - defect cases."""

    defectFlux = 100000  # Flux to use for simulated defect.

    def setUp(self):
        """Generate a mock exposure/camera to test."""
        config = isrMock.IsrMockConfig()
        config.isTrimmed = True
        config.doGenerateImage = True
        config.doAddFringe = False
        config.doAddSource = False
        config.doAddSky = True
        config.doAddOverscan = False
        config.doAddCrosstalk = False
        config.doAddBias = False
        config.doAddDark = False
        config.doAddFlat = False
        config.doAddFringe = False

        config.skyLevel = 1000
        config.rngSeed = 12345
        self.inputExp = isrMock.IsrMock(config=config).run()

        # These are simulated defects
        self.inputExp.getImage().getArray()[0, 0] = -1.0 * self.defectFlux
        self.inputExp.getImage().getArray()[40, 50] = self.defectFlux
        self.inputExp.getImage().getArray()[75, 50] = np.nan

        updateMockExp(self.inputExp, addCR=False)

        self.inputExp.getMask().getArray()[0, 0] = 1
        self.inputExp.getMask().getArray()[40, 50] = 1
        self.inputExp.getMask().getArray()[75, 50] = 1

        self.camera = isrMock.IsrMock().getCamera()

    def test_defects(self):
        """Test a subset of the output values to identify that the
        image stat methods haven't changed.
        """
        config = cpVerify.CpVerifyDefectsConfig()
        config.numSigmaClip = 3.0
        task = cpVerify.CpVerifyDefectsTask(config=config)
        results = task.run(self.inputExp, self.camera)
        defectStats = results.outputStats

        self.assertEqual(defectStats['AMP']['C:0,0']['DEFECT_PIXELS'], 53)
        self.assertEqual(defectStats['AMP']['C:0,0']['OUTLIERS'], 17)
        self.assertEqual(defectStats['AMP']['C:0,0']['STAT_OUTLIERS'], 3)
        self.assertAlmostEqual(defectStats['AMP']['C:0,0']['MEDIAN'], 999.466, 4)
        self.assertAlmostEqual(defectStats['AMP']['C:0,0']['STDEV'], 30.96303, 4)
        self.assertAlmostEqual(defectStats['AMP']['C:0,0']['MIN'], 881.56146, 4)
        self.assertAlmostEqual(defectStats['AMP']['C:0,0']['MAX'], 1124.19934, 4)

        self.assertEqual(defectStats['AMP']['C:0,0']['UNMASKED_MIN'], -1.0 * self.defectFlux, 4)
        self.assertEqual(defectStats['AMP']['C:0,0']['UNMASKED_MAX'], self.defectFlux, 4)

        self.assertIn(defectStats['SUCCESS'], [True, False])


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
