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
import lsst.cp.verify as cpVerify

from lsst.ip.isr import CrosstalkCalib


class VerifyCrosstalkTestCase(lsst.utils.tests.TestCase):
    """Unit test for stats code - crosstalk cases.
    """

    def test_crosstalk(self):
        """Test a subset of the output values to identify that the
        image stat methods haven't changed.
        """
        crosstalk = CrosstalkCalib(nAmp=2)
        crosstalk.coeffs = np.array([[0.0, 1e-4], [1e-4, 0.0]])
        crosstalk.coeffErr = np.array([[0.0, 1e-5], [1e-5, 0.0]])
        crosstalk.coeffNum = np.array([[0, 100], [100, 0]])
        crosstalk.coeffValid = np.array([[False, True], [True, False]], dtype=bool)

        config = cpVerify.CpVerifyCrosstalkConfig()
        task = cpVerify.CpVerifyCrosstalkTask(config=config)
        results = task.run(crosstalk)
        crosstalkStats = results.outputStats

        self.assertEqual(crosstalkStats['DETECTOR']['N_AMP'], crosstalk.nAmp)
        self.assertEqual(crosstalkStats['DETECTOR']['N_VALID'], 2)

        self.assertFalse(crosstalkStats['SUCCESS'])
        self.assertFalse(crosstalkStats['VERIFY']['NO_SIGNIFICANT_DETECTION'])

    def test_crosstalkNull(self):
        """Test a subset of the output values to identify that the
        image stat methods haven't changed.
        """
        crosstalk = CrosstalkCalib(nAmp=2)
        crosstalk.coeffs = np.array([[0.0, 1e-6], [1e-6, 0.0]]),
        crosstalk.coeffErr = np.array([[0.0, 1e-5], [1e-5, 0.0]]),
        crosstalk.coeffNum = np.array([[0, 100], [100, 0]]),
        crosstalk.coeffValid = np.array([[False, False], [False, False]], dtype=bool)

        config = cpVerify.CpVerifyCrosstalkConfig()
        task = cpVerify.CpVerifyCrosstalkTask(config=config)
        results = task.run(crosstalk)
        crosstalkStats = results.outputStats

        self.assertEqual(crosstalkStats['DETECTOR']['N_AMP'], crosstalk.nAmp)
        self.assertEqual(crosstalkStats['DETECTOR']['N_VALID'], 0)

        self.assertTrue(crosstalkStats['SUCCESS'])
        self.assertTrue(crosstalkStats['VERIFY']['NO_SIGNIFICANT_DETECTION'])


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
