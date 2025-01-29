BPS Template Files for Calibration Construction
===============================================

This directory contains BPS template files for calibration verification.  These
are meant to conform to the naming conventions laid out in
[DMTN-222](https://dmtn-222.lsst.io).

Overview
--------

These template files are only part of calibration verification. They are meant
to be flexible enough to use for both test verification into a user collection,
as well as final verification.

The templates rely on a lot of exported environment variables, defined
below. Note that because the bps file does not share an env with the caller all
the env vars must be explicitly exported.

More information on bps is [here](https://pipelines.lsst.io/modules/lsst.ctrl.bps/quickstart.html).

Environment Variables
---------------------

These are the environment variables that must be used. Examples are substituted in.

* `export USER_CALIB_PREFIX=""` or `export USER_CALIB_PREFIX=u/erykoff/`: Set to an empty string for official calibrations, or the user collection prefix (**with trailing slash**).
* `export INSTRUMENT=LSSTComCam`: The name of the instrument.
* `export TICKET=DM-46562`: The name of the ticket assocated with the calib construction.
* `export REPO=/repo/main`: The name of the butler repository to generate calibs.
* `export RAW_COLLECTION=LSSTComCam/raw/all`: The name of the raw data collection.
* `export CALIB_COLLECTIONS=LSSTComCam/calib/DM-46825`: Comma-separated list of curated or previously generated calibration collections to use as a starting point.
* `export TAG=newCalibs`: A human-readable tag to help indicate why a set of calibs were built (should also be findable from the ticket name).
* `export VERIFY_RERUN=20250122a`: The rerun name to indicate when the verification was run.
* `export SELECTION_BIAS_V="instrument='LSSTComCam' and selection_string"`: The selection of raws to verify the bias.
* `export SELECTION_DARK_V="instrument='LSSTComCam' and selection_string"`: The selection of raws to verify the dark.
* `export SELECTION_PTC_V="instrument='LSSTComCam' and selection_string"`: The selection of raws to verify the PTC.
* `export SELECTION_PTC_LINEARIZER_V=$SELECTION_PTC_V`: The selection of raws to verify the linearizer; usually will be the same as the PTC selection.
* `export SELECTION_PTC_BFK_V=$SELECTION_PTC_V`: The selection of raws to verify the brighter-fatter kernel; usually will be the same as the PTC selection.
* `export SELECTION_FLAT_g_V="instrument='LSSTComCam' and selection_string": The selection of raws to verify the g-band flat. See below for additional info.

Checking Environment Variables
------------------------------

Before submitting a bps file, it can be helpful to check that you have all the
correct environment variables set.  This can be mostly accomplished with the
following bash command (using PTC as an example):

```
cat $CP_VERIFY_DIR/bps/templates/bps_verify_ptc.yaml | envsubst
```

This will render all the env vars as bash would render them. Unfortunately, bps
uses python `os.path.expandvars()` which is not as sophisticated. The most
relevant limitation is that `envsubst` will replace an unset env var with a
empty string, while `expandvars()` will leave that raw
(e.g. `${USER_CALIB_PREFIX}` can be left verbatim). Nevertheless, this is a
useful sanity check before submission.

A Note on Flat Verification
---------------------------

Due to limitations in BPS environment variable expansion, we have one template
file for each band. There are template files for each of ugrizy. Note that it is
**strongly** recommended that the relevant physical filter be used in the
`selection_string` constraint on each flat selection.
