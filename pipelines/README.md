# Pipeline Definitions

This directory contains pipeline definition YAML files which are used to verify calibration products with the LSST Science Pipelines.
The pipelines defined here are science ready and come in two flavors: generic ("ingredients") and camera-specific (within appropriately named sub-directories).
Use of camera-specific pipelines is strongly encouraged where possible as they are optimized for the particular characteristics of that camera.

The pipelines defined here tend to import other pipelines, including ingredient pipelines in the [_ingredients](_ingredients) directory.
To expand a pipeline YAML and resolve such imports for the purposes of visualizing it, the `pipetask build` command can be used.
For example, to visualize the verification of a bias for the [LATISS camera pipeline](https://github.com/lsst/cp_verify/blob/main/pipelines/Latiss/VerifyBias.yaml) pipeline, run:

```bash
pipetask build \
-p $CP_VERIFY_DIR/pipelines/Latiss/VerifyBias.yaml \
--show pipeline
```
