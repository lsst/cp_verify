.. py:currentmodule:: lsst.cp.verify

#####################################
Verifying calibrations with cp_verify
#####################################

`cp_verify` is designed to make the verification of calibrations against the criteria listed in `DMTN-101 <https://dmtn-101.lsst.io>`_ as simple as possible.  However, this has resulted in the tasks being rather abstract, with no clear connection between steps.  This document should help explain how the tests are run, how those results are stored, and how they are passed to `analysis_tools` for conversion into metrics and plots.

Last updated: 2024-04-19.

################
Basic principles
################

The `cp_verify` process relies on `pipeline <_pipe_base_creating_pipeline>_` definitions linking a series of `pipetask <_pipe_base-creating-a-pipelinetask>_` commands to process a number of input exposures using the calibration to be verified.  These processed exposures then have quantities measured by the appropriate `cp_verify` task (generally a subclass of the `CpVerifyStatsTask`), writing out results containing those quantities, and whether the measurements are within the specifications listed in DMTN-101.  These tasks operate on data in a per-detector and per-exposure manner, so to obtain a determination if the calibration has met the requirements, they are first accumulated into a per-exposure dataset that summarizes the results over all detectors using one of the `CpVerifyExpMergeTask` subclasses.  Next, these per-exposure datasets are again accumulated by another such task into a per-run dataset that merges all exposures together.  This final dataset should summarize the entirety of the verification process, and provide the information needed for calibration evaluation and acceptance.

This structure is shared by most of the verification tasks, however some calibrations, such as the photon transfer curve (PTC), can be studied without application to input exposures.  These calibrations are directly studied by a subclass of the `CpVerifyCalibTask`.  These subclasses again check quantities (this time taken directly from the calibration) against the specifications in DMTN-101, and accumulate per-exposure and per-run datasets.  The two base classes are similar in structure, however they differ due to the need for different <lsst.daf.butler.Connection> definitions.

The final accumulated results from either style of pipeline can then be passed to the appropriate `analysis_tools` tasks to create the metrics and plots distributed to their final destination.  This section will be updated when this `analysis_tools` step has been fully determined.

########################
A note about collections
########################

The majority of calibrations being checked by `cp_verify` have not been accepted and certified for a specific time range.  To ensure that the correct calibration is being verified, the final output :ref:`RUN collection <lsst.daf.butler.CollectionType.RUN>` should be the first input collection passed to the pipeline.  Failing to include this collection (or having that collection come after other collections which may contain a combined calibration of the same kind) will result in the verification processing either failing or silently proceeding with a different combined calibration.

###############################
Example: verification of a bias
###############################

To investigate the structure of these tasks, their configurations, and the expected results, let us consider the operations needed for the verification of a combined bias.  Following the pipeline definition, the first step of bias verification is the processing of all input exposures through `IsrTask` using the configuration overrides present in the pipeline definition.  These overrides ensure the tasks applies all stages of instrument signal removal up to and including the bias subtraction, but none of the image processing steps that follow this step.  As we wish to evaluate the quality of the combined bias calibration, the input exposures should be bias frames (zero second exposures taken with the shutter closed).  Applying the combined bias to the bias frames should remove all signal from the image, so the expected output is an image with zero mean and noise consistent with the read noise of the detector (or, in most cases, the read noise of the detector segment corresponding to the read-out amplifiers).  As the images are already in memory at this stage, the `lsst.ip.isr.IsrStatisticsTask` is used to do measurements that need the calibration data, the overscan regions, or some other product that is more easily available through the `IsrTask` connections.

The outputs of `IsrTask` are best thought of as "residual" images showing what remaining signal is still present.  As a different configuration is used for these processings than would be in a standard science analysis, these residual images use a different dataset type (``verifyBiasProc`` for bias verification) to distinguish them from the standard ``postISRCCD` products.

These residual images are the main input to `CpVerifyBiasTask`, a subclass of `CpVerifyStatsTask`, along with the `IsrTask` task metadata (for historical reasons; this was used as the source of the empirical read noise estimate), and any ``isrStatistics`` results.  The parent class has a handful of other input connections, allowing catalog results to be verified, and uncorrected versions of the input exposure and catalog (which have not had the calibration applied).  These connections are discussed below as they apply to particular calibration products.

For simple measurements, the `CpVerifyBiasConfig` defines a set of dictionary configuration entries that have string keys that correspond to a measurement name that will appear in the output results, with a value that maps in most cases maps to an `lsst.afw.math.statistics.Property` indicating what mathematical operation will be used to perform that measurement.  If these dictionary configurations have non-zero length, then they are used by the task to determine which methods will be used to do the measurement.  These configuration options include:

* ``imageStatKeywords``:  these options define the tests that will be run on each amplifier segment in the image being processed.  All following ``ImageStatKeyword`` entries follow this same behavior.

* ``unmaskedImageStatKeywords``:  as above, these options define the tests run ignoring all mask planes present.

* ``uncorrectedImageStatKeywords``:  these options define the tests that will be run using the ``uncorrectedExp`` input connection.

* ``crImageStatKeywords``:  these options define the tests that will be run after cosmic ray rejection has been applied.

* ``normImageStatKeywords``:  these options define the tests that will be run after normalizing the image by the exposure time.

* ``detectorStatKeywords``:  as all previous image analysis is done on a per-amplifier basis, these options define tests to be run on the full detector.  As these will nearly always depend on the calibration under consideration, there is no defined method in the base class, so subclasses will need to implement it.

* ``metadataStatKeywords``:  these options define the tests that will be run on the input task_metadata input.  These tests need to be implemented by subclasses.

* ``catalogStatKeywords``:  these options define the tests that will be run on the input catalog data.  As before, these tests need to be implemented by subclasses.

The final two configuration options that control the behavior of the verification task are the ``useIsrStatistics`` configuration option, which defines if there are results from the ``isrStatistics`` input that should be considered, and the ``hasMatrixCatalog`` option, which indicates if a matrix catalog output will be constructed.

The output results are stored in one of three formats.  The ``outputStats`` (``verifyBiasDetStats`` for the bias example) contains the legacy output as a YAML output of nested dictionaries.  The ``outputResults`` (``verifyBiasDetResults``) is the newer flat table catalog containing the same results as the ``outputStats``, but presented with each test result and other quantities in each column, and each row containing results for either an amplifier segment or the full detector.  The ``outputMatrix`` (``verifyBiasDetMatrix``) contains a catalog containing results that can be used to construct one or many matrices of values, with each row of the catalog containing the values for a single element of the final matrices.  These dataset types have been chosen to indicate that these data products contain the results from a single detector in the exposure.

The main task checks for config options that have non-zero length, and call the appropriate method to generate the results indicated in those configuration options.  Once all measurements are done, the combined ``outputStats`` dictionary data is passed to the ``verify`` method.  This method must always be implemented by the subclass, as it contains the code to perform the checks described in DMTN-101 to evaluate the calibration quality.  Generally, each test name from the configuration options are evaluated and stored as an entry in a "verify" dictionary with the same key.  As much as is possible, the test checks in the ``verify`` method should indicate which test from DMTN-101 is being checked, and return a simple boolean.

During development, a special ``FORCE_FAILURE`` statistic has been used to indicate that some operation (generally the cosmic ray rejection) has failed, and that despite the results of the other tests, this analysis should be considered a failure.  As cosmic ray rejection has been improved, this is no longer a common result, but may be useful for understanding existing analysis.

After evaluating the DMTN-101 test success, the nested dictionaries that are returned in the ``outputStats`` are passed to the ``repackStats`` method to construct the catalog tables for the ``outputResults`` and ``outputMatrix``.  In the future, we will likely deprecate the ``outputStats`` product entirely, and refactor the tasks to directly write to these tables.  As the data stored in the ``outputResults`` product is what will be passed to `analysis_tools`, it must contain the column names expected by those tasks.

The merging process involves the ``Stats``, ``Results``, and ``Matrix`` products from multiple detectors being merged into a single output per exposure (e.g. ``verifyBiasExpResults``) and those per exposure products merged into the final result summarizing the full processing run (``verifyBiasRunResults``).  The ``Stats`` product is sparse, only containing information about failures and any new results calculated on a per-exposure level.  In contrast to this, the ``Results`` and ``Matrix`` are concatenations of all the inputs, such that the final ``verifyBiasRunResults`` catalog contains measurements across all amplifiers from all detectors from all exposures processed.

As the merging process generally requires fewer per-calibration type special cases, the task is intended to be as generic as possible, using the task configuration fields to control the merging process.  The important configuration parameters are:

* ``statKeywords``:  as was done in the verify task, this is a dictionary with test names as keys, and values indicating the test to be performed.  Using this requires a subclass to implement the ``calcStatistics`` method.

* ``hasMatrixCatalog``, ``hasInputResults``:  these options indicate whether these input connections will be populated, and if they are populated, that they will need to be merged.  With these options set to false, only the legacy ``Stats`` products will be used.

* ``mergeDimension``:  the same base class is used for merging per-detector results and for merging per-exposure results.  This config option controls which dimension is being summed by the associated task.

* ``stageName``:  the stage name must be set here, to ensure that log messages and catalog columns correctly indicate the product being verified.

############################
Calibration specific details
############################


