#!/usr/bin/env python

#
# LSST Data Management System
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
"""Test cases for cp_verify pipelines."""

import unittest

from lsst.pipe.base import Pipeline, PipelineGraph
from lsst.resources import ResourcePath
import lsst.utils

try:
    import lsst.obs.lsst
    has_obs_lsst = True
except ImportError:
    has_obs_lsst = False

try:
    import lsst.obs.subaru
    has_obs_subaru = True
except ImportError:
    has_obs_subaru = False

try:
    import lsst.obs.decam
    has_obs_decam = True
except ImportError:
    has_obs_decam = False

PIPELINE_URI = ResourcePath("eups://cp_verify/pipelines/", forceDirectory=True)


class VerifyPipelinesTestCase(lsst.utils.tests.TestCase):
    """Test case for building the pipelines."""

    def _get_pipelines(self, exclude=[]):
        pipelines = {
            "verifyBfk.yaml",
            "verifyBias.yaml",
            "verifyCrosstalk.yaml",
            "verifyDark.yaml",
            "verifyDefectsIndividual.yaml",
            "verifyDefects.yaml",
            "verifyFlat.yaml",
            "verifyScience.yaml",
            # Old pipeline name.
            "verifyGain.yaml",
            # New pipeline name.
            "verifyGainFromFlatPairs.yaml",
            "verifyLinearizer.yaml",
            "verifyPtc.yaml",
            "verifyIlluminationCorrection.yaml",
        }

        for ex in exclude:
            pipelines.remove(ex)

        return pipelines

    def _check_pipeline(self, pipeline_file: ResourcePath):
        # Confirm that the file is there.
        self.assertTrue(pipeline_file.exists(), msg=f"Could not find {pipeline_file}")

        # The following loads the pipeline and confirms that it can parse all
        # the configs.
        try:
            pipeline = Pipeline.from_uri(pipeline_file)
            graph = pipeline.to_graph()
        except Exception as e:
            raise RuntimeError(f"Could not process {pipeline_file}") from e

        self.assertIsInstance(graph, PipelineGraph)

    def test_ingredients(self):
        """Check that all pipelines in pipelines/_ingredients are tested."""
        ingredient_files = ResourcePath.findFileResources(
            [PIPELINE_URI.join("_ingredients")], file_filter=r".*\.yaml$"
        )
        # The *LSST.yaml pipelines are imported by LATISS/LSSTComCam/LSSTCam
        # and are not tested on their own.
        ingredients = set(
            [pipeline.basename() for pipeline in ingredient_files if "LSST.yaml" not in pipeline.path]
        )
        # The _ingredients/verifyGainFromFlatPairs.yaml becomes
        # verifyGain.yaml in older pipelines for compatibility.
        expected = self._get_pipelines()
        # The _ingredients/verifyGainFromFlatPairs.yaml becomes
        # verifyGain.yaml in older pipelines for compatibility.
        expected.remove("verifyGain.yaml")
        # There is no regular verifyScience.yaml (non-LSST)
        expected.remove("verifyScience.yaml")
        # The Illumination Correction pipeline is only an "LSST" version.
        expected.remove("verifyIlluminationCorrection.yaml")
        self.assertEqual(ingredients, expected)

    def test_cameras(self):
        """Check that all the cameras in pipelines are tested."""
        _, paths, _ = next(PIPELINE_URI.walk())
        expected = {
            "DECam",
            "HSC",
            "_ingredients",
            "LATISS",
            "LSSTCam",
            "LSSTCam-imSim",
            "LSSTComCam",
            "LSSTComCamSim",
        }
        self.assertEqual(set(paths), expected)

    @unittest.skipIf(not has_obs_lsst, reason="Cannot test LATISS pipelines without obs_lsst")
    def test_latiss_pipelines(self):
        latiss_uri = PIPELINE_URI.join("LATISS", forceDirectory=True)
        for pipeline in self._get_pipelines(exclude=[
                # The old pipeline name should be excluded.
                "verifyGain.yaml",
                # The following tasks are not part of the new pipelines.
                "verifyDefectsIndividual.yaml",
                # The following tasks will be added in the future.
                "verifyCrosstalk.yaml",
                "verifyBfk.yaml",
                "verifyIlluminationCorrection.yaml",

        ]):
            self._check_pipeline(latiss_uri.join(pipeline))

    @unittest.skipIf(not has_obs_lsst, reason="Cannot test LSSTCam pipelines without obs_lsst")
    def test_lsstcam_pipelines(self):
        lsstcam_uri = PIPELINE_URI.join("LSSTCam", forceDirectory=True)
        for pipeline in self._get_pipelines(
                exclude=[
                    # These are renamed/not used in the new pipelines.
                    "verifyGain.yaml",
                    "verifyDefectsIndividual.yaml",
                    # These are not used yet.
                    "verifyCrosstalk.yaml",
                    "verifyIlluminationCorrection.yaml",
                ]):
            self._check_pipeline(lsstcam_uri.join(pipeline))

    @unittest.skipIf(not has_obs_lsst, reason="Cannot test LSSTCam-imSim pipelines without obs_lsst")
    def test_lsstcam_imsim_pipelines(self):
        sim_uri = PIPELINE_URI.join("LSSTCam-imSim", forceDirectory=True)
        for pipeline in self._get_pipelines(
            exclude=[
                "verifyGainFromFlatPairs.yaml",
                "verifyScience.yaml",
                "verifyIlluminationCorrection.yaml",
                "verifyCrosstalk.yaml",
            ],
        ):
            self._check_pipeline(sim_uri.join(pipeline))

    @unittest.skipIf(not has_obs_lsst, reason="Cannot test LSSTComCam pipelines without obs_lsst")
    def test_lsstcomcam_pipelines(self):
        comcam_uri = PIPELINE_URI.join("LSSTComCam", forceDirectory=True)
        for pipeline in self._get_pipelines(
            exclude=[
                # These are renamed/not used in the new pipelines.
                "verifyGain.yaml",
                "verifyDefectsIndividual.yaml",
                # These are not used yet.
                "verifyCrosstalk.yaml",
            ],
        ):
            self._check_pipeline(comcam_uri.join(pipeline))

    @unittest.skipIf(not has_obs_lsst, reason="Cannot test LSSTComCamSim pipelines without obs_lsst")
    def test_lsstcomcamsim_pipelines(self):
        comcam_sim_uri = PIPELINE_URI.join("LSSTComCamSim", forceDirectory=True)
        for pipeline in self._get_pipelines(
            exclude=[
                # These are renamed/not used in the new pipelines.
                "verifyGain.yaml",
                "verifyDefectsIndividual.yaml",
                # These are not valid for LSSTComCamSim.
                "verifyCrosstalk.yaml",
                "verifyLinearizer.yaml",
                "verifyIlluminationCorrection.yaml",

            ],
        ):
            self._check_pipeline(comcam_sim_uri.join(pipeline))

    @unittest.skipIf(not has_obs_decam, reason="Cannot test DECam pipelines without obs_decam")
    def test_decam_pipelines(self):
        decam_uri = PIPELINE_URI.join("DECam", forceDirectory=True)
        for pipeline in self._get_pipelines(
            exclude=[
                "verifyGainFromFlatPairs.yaml",
                "verifyScience.yaml",
                "verifyIlluminationCorrection.yaml",
                "verifyCrosstalk.yaml",

            ],
        ):
            self._check_pipeline(decam_uri.join(pipeline))

    @unittest.skipIf(not has_obs_subaru, reason="Cannot test HSC pipelines without obs_subaru")
    def test_hsc_pipelines(self):
        hsc_uri = PIPELINE_URI.join("HSC", forceDirectory=True)
        for pipeline in self._get_pipelines(
            exclude=[
                "verifyGainFromFlatPairs.yaml",
                "verifyScience.yaml",
                "verifyIlluminationCorrection.yaml",
                "verifyCrosstalk.yaml",

            ],
        ):
            self._check_pipeline(hsc_uri.join(pipeline))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
