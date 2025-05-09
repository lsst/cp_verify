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

from astropy.table import Table

import lsst.utils.tests
import lsst.afw.cameraGeom as cameraGeom
import lsst.cp.verify as cpVerify


class MergeResultsTestCase(lsst.utils.tests.TestCase):
    """Unit test for merge code.
    """

    def setUp(self):
        # Generate a camera.
        self.camera = cameraGeom.testUtils.CameraWrapper(isLsstLike=False).camera

        # Fix the seed so we can expect consistent results.
        np.random.seed(12345)

        # Make first exposure
        self.detectorResults1 = list()
        self.detectorDimensions1 = list()
        for detector in self.camera:
            dimensions = {'exposure': 1, 'detector': detector.getId(), 'instrument': 'testCam'}
            self.detectorDimensions1.append(dimensions)
            self.detectorResults1.append(self.generateDetectorResults(detector))

        # Make second exposure
        self.detectorResults2 = list()
        self.detectorDimensions2 = list()
        for detector in self.camera:
            dimensions = {'exposure': 2, 'detector': detector.getId(), 'instrument': 'testCam'}
            self.detectorDimensions2.append(dimensions)
            self.detectorResults2.append(self.generateDetectorResults(detector, forceFailure=True))

    def test_merging(self):
        """Generate simulated verify statistics, and prove they merge
        correctly.
        """
        # Generate table versions of the detector results.
        verifyStatsTask = cpVerify.CpVerifyStatsTask()
        expTables1 = [Table(verifyStatsTask.repackStats(rr, dd)[0]) for rr, dd in
                      zip(self.detectorResults1, self.detectorDimensions1)]
        expTables2 = [Table(verifyStatsTask.repackStats(rr, dd)[0]) for rr, dd in
                      zip(self.detectorResults2, self.detectorDimensions2)]

        # Ensure we're merging the tables as well.
        expMergeConfig = cpVerify.CpVerifyExpMergeTask.ConfigClass()
        expMergeConfig.hasInputResults = True

        # Do the merge.
        expMergeTask = cpVerify.CpVerifyExpMergeTask(config=expMergeConfig)
        expResults1 = expMergeTask.run(self.detectorResults1, self.detectorDimensions1,
                                       self.camera, inputResults=expTables1)
        expResults2 = expMergeTask.run(self.detectorResults2, self.detectorDimensions2,
                                       self.camera, inputResults=expTables2)
        self.assertTrue(expResults1.outputStats['SUCCESS'])
        self.assertFalse(expResults2.outputStats['SUCCESS'])  # Expected to have one failure
        self.assertEqual(expResults2.outputStats['R:1,0 S:1,0']['FAILURES'][0], "C:0,2 SIGMA_TEST")

        # Prep inputs to full calibration run merge.
        inputStats = [expResults1.outputStats, expResults2.outputStats]
        inputDims = [{'exposure': 1}, {'exposure': 2}]

        runMergeTask = cpVerify.CpVerifyRunMergeTask()

        runResults = runMergeTask.run(inputStats, inputDims, self.camera)
        self.assertFalse(runResults.outputStats['SUCCESS'])

        # We know this one failed.
        failureList = runResults.outputStats[2]['FAILURES']
        self.assertEqual(len(failureList), 1)
        self.assertEqual(failureList[0], "R:1,0 S:1,0 C:0,2 SIGMA_TEST")

    def generateDetectorResults(self, detector, mean=10.0, sigma=1.2, forceFailure=False):
        """Make the simulated verify statistics.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.Detector`
            Detector geometry to make statistics for.
        mean : `float`
            Center of the random normal distribution to pull
            statistics from.
        sigma : `float`
            Sigma of the normal distribution.
        forceFailure : `bool`
            If True, force the "verify" to fail.

        Returns
        -------
        detectorStats : `dict` [`str`, `dict`]
            Nested dictionary of verify statistics.
        """

        valueNames = ['MEAN', 'SIGMA', 'VALUE']

        ampDict = dict()
        detDict = dict()
        catDict = dict()
        success = True
        for amp in detector.getAmplifiers():
            ampName = amp.getName()
            ampDict[ampName] = dict()
            ampDict[ampName]['VECTOR'] = []
            for value in valueNames:
                ampDict[ampName][value] = np.random.normal(10.0, 1.2)
            ampDict[ampName]['VECTOR'].append(np.random.uniform(size=len(detector)))

        ampVerify = dict()
        detVerify = dict()
        catVerify = dict()
        expVerify = dict()
        for ampName in ampDict:
            ampVerify[ampName] = dict()
            ampSuccess = True
            for value in valueNames:
                if forceFailure and (value == "SIGMA"
                                     and ampName == "C:0,2"
                                     and detector.getName() == 'R:1,0 S:1,0'):
                    ampVerify[ampName][value + "_TEST"] = False
                    ampSuccess = False
                    success = False
                else:
                    ampVerify[ampName][value + "_TEST"] = True
            ampVerify[ampName]['SUCCESS'] = ampSuccess

        return {'SUCCESS': success, 'AMP': ampDict, 'DET': detDict, 'CATALOG': catDict,
                'VERIFY': {'AMP': ampVerify, 'DET': detVerify, 'CATALOG': catVerify, 'EXP': expVerify}}


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
