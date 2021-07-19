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

import pprint
import defaultdict
import numpy as np

def interactiveBlock(description, detStats):
    """Handle interaction during loops.

    Parameters
    ----------
    description : `str`
        Description to use in prompt.
    detStats : `dict` [`str` `dict`]
        Nested dictionaries containing the current set of statistics.

    Returns
    -------
    continueDisplay : `bool`
        Whether further examples should continue.
    nSkip : `int`
        Number of further examples to skip.
    """
    while True:
        ans = input(f"{description} Continue? [c, q, p, #]").lower()
        if ans in ("", "c"):
            break
        elif ans in ("q", "x"):
            return False, 0
        elif ans in ('p', ):
            pprint(detStats)
            break
        elif ans.isnumeric():
            return True, int(ans)
    return True, 0

# Taken from faro
def calcQuartileClippedStats(dataArray, nSigmaToClip=3.0):
    """Calculate the quartile-based clipped statistics of a data array.
    The difference between quartiles[2] and quartiles[0] is the interquartile
    distance.  0.74*interquartileDistance is an estimate of standard deviation
    so, in the case that ``dataArray`` has an approximately Gaussian
    distribution, this is equivalent to nSigma clipping.
    Parameters
    ----------
    dataArray : `list` or `numpy.ndarray` of `float`
        List or array containing the values for which the quartile-based
        clipped statistics are to be calculated.
    nSigmaToClip : `float`, optional
        Number of \"sigma\" outside of which to clip data when computing the
        statistics.
    Returns
    -------
    result : `lsst.pipe.base.Struct`
        The quartile-based clipped statistics with ``nSigmaToClip`` clipping.
        Atributes are:
        ``median``
            The median of the full ``dataArray`` (`float`).
        ``mean``
            The quartile-based clipped mean (`float`).
        ``stdDev``
            The quartile-based clipped standard deviation (`float`).
        ``rms``
            The quartile-based clipped root-mean-squared (`float`).
        ``clipValue``
            The value outside of which to clip the data before computing the
            statistics (`float`).
        ``goodArray``
            A boolean array indicating which data points in ``dataArray`` were
            used in the calculation of the statistics, where `False` indicates
            a clipped datapoint (`numpy.ndarray` of `bool`).
    """
    quartiles = np.percentile(dataArray, [25, 50, 75])
    assert len(quartiles) == 3
    median = quartiles[1]
    interQuartileDistance = quartiles[2] - quartiles[0]
    clipValue = nSigmaToClip*0.74*interQuartileDistance
    good = np.logical_not(np.abs(dataArray - median) > clipValue)
    quartileClippedMean = dataArray[good].mean()
    quartileClippedStdDev = dataArray[good].std()
    quartileClippedRms = np.sqrt(np.mean(dataArray[good]**2))

    return Struct(
        median=median,
        mean=quartileClippedMean,
        stdDev=quartileClippedStdDev,
        rms=quartileClippedRms,
        clipValue=clipValue,
        goodArray=good,
    )

# Taken from pipe_analysis
def plotCameraOutline(axes, camera, ccdList, color="k", fontSize=6, metricPerCcdDict=None,
                      metricStr="", fig=None, metricSigmaRange=4.0):
    """Plot the outline of the camera ccds highlighting those with data.

    Parameters
    ----------
    axes : `matplotlib.axes._axes.Axes`
        Particular matplotlib axes on which to plot the tract outline.
    camera : `lsst.afw.cameraGeom.Camera`, optional
        The ``camera`` associated with the dataset (used to label the plot with
        the camera's name).
    ccdList : `list` of `int`, optional
        List of ccd IDs with data to be plotted.
    fontSize : `int`, optional
        Font size for plot labels.
    metricPerCcdDict : `dict` of `float`, optional
        Dictionary of per patch metric averages; {ccdId: metricValue}.  If
        provided, these values will be used to color-code the camera outline
        plot.
    metricStr : `str`, optional
        String representing the computed metric values provided in
        ``metricPerCcdDict``.
    fig : `matplotlib.figure.Figure`, optional
        The figure on which to add the per-ccd metric info (required to add
        the colorbar).
    metricSigmaRange : `float`, optional
        Number of sigma to make the +/- range for the metric colorbar.
    """
    if metricPerCcdDict:
        if fig is None:
            raise RuntimeError("Must supply the matplotlib fig if color-coding by metric-per-ccd")
    axes.tick_params(which="both", direction="in", labelleft=False, labelbottom=False)
    axes.locator_params(nbins=6)
    axes.ticklabel_format(useOffset=False)
    camRadius = max(camera.getFpBBox().getWidth(), camera.getFpBBox().getHeight())/2
    camRadius = np.round(camRadius, -1)
    camLimits = np.round(1.25*camRadius, -1)
    intCcdList = [int(ccd) for ccd in ccdList]
    prop_cycle = plt.rcParams["axes.prop_cycle"]
    colors = prop_cycle.by_key()["color"]
    colors.pop(colors.index("#7f7f7f"))  # get rid of the gray one as is doesn't contrast well with white
    colors.append("gold")
    hasRotatedCcds = False
    for ccd in camera:
        if ccd.getOrientation().getNQuarter() != 0:
            hasRotatedCcds = True
            break
    if metricPerCcdDict:  # color-code the ccds by the per-ccd metric measurement
        cmap = plt.cm.viridis
        metricPerCcdArray = np.fromiter(metricPerCcdDict.values(), dtype="float32")
        finiteMetrics = np.isfinite(metricPerCcdArray)
        clippedStats = calcQuartileClippedStats(metricPerCcdArray[finiteMetrics], nSigmaToClip=5.0)
        vMin = clippedStats.mean - metricSigmaRange*clippedStats.stdDev
        vMax = clippedStats.mean + metricSigmaRange*clippedStats.stdDev
        vMax = max(abs(vMin), vMax) if vMax > 0 else vMax  # Make range symmetric about 0 if it crosses 0
        vMin = -vMax if vMax > 0 else vMin
        cmapBins = np.linspace(vMin, vMax, cmap.N - 1)
    for ic, ccd in enumerate(camera):
        ccdCorners = ccd.getCorners(cameraGeom.FOCAL_PLANE)
        if ccd.getType() == cameraGeom.DetectorType.SCIENCE:
            axes.add_patch(patches.Rectangle(ccdCorners[0], *list(ccdCorners[2] - ccdCorners[0]),
                           facecolor="none", edgecolor="k", ls="solid", lw=0.5, alpha=0.5))
        if ccd.getId() in intCcdList:
            if metricPerCcdDict is None:
                if hasRotatedCcds:
                    nQuarter = ccd.getOrientation().getNQuarter()
                    fillColor = colors[nQuarter%len(colors)]
                elif ccd.getName()[0] == "R":
                    try:
                        fillColor = colors[(int(ccd.getName()[1]) + int(ccd.getName()[2]))%len(colors)]
                    except Exception:
                        fillColor = colors[ic%len(colors)]
                else:
                    fillColor = colors[ic%len(colors)]
            else:
                cmapBinIndex = np.digitize(metricPerCcdDict[str(ccd.getId())], cmapBins)
                fillColor = cmap.colors[cmapBinIndex]
            ccdCorners = ccd.getCorners(cameraGeom.FOCAL_PLANE)
            axes.add_patch(patches.Rectangle(ccdCorners[0], *list(ccdCorners[2] - ccdCorners[0]),
                                             fill=True, facecolor=fillColor, edgecolor="k",
                                             ls="solid", lw=1.0, alpha=0.7))
    axes.set_xlim(-camLimits, camLimits)
    axes.set_ylim(-camLimits, camLimits)
    if camera.getName() == "HSC":
        for x, y, t in ([-1, 0, "N"], [0, 1, "W"], [1, 0, "S"], [0, -1, "E"]):
            axes.text(1.085*camRadius*x, 1.085*camRadius*y, t, ha="center", va="center",
                      fontsize=fontSize - 1)
            axes.text(-0.82*camRadius, 1.04*camRadius, "%s" % camera.getName(), ha="center",
                      fontsize=fontSize, color=color)
    else:
        axes.text(0.0, 1.04*camRadius, "%s" % camera.getName(), ha="center", fontsize=fontSize,
                  color=color)
    if metricPerCcdDict:
        mappable = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vMin, vmax=vMax))
        mappable._A = []  # fake up the array of the scalar mappable. Urgh...
        axesBbox = axes.get_position()
        caxDim = [axesBbox.xmin, 0.95*axesBbox.ymin, axesBbox.width, 0.07*axesBbox.height]
        cax = fig.add_axes(caxDim)
        cax.tick_params(labelsize=fontSize - 1)
        cbar = fig.colorbar(mappable, cax=cax, orientation="horizontal")
        cbar.set_label(label=metricStr, size=fontSize)


def plotFailures(runStats, camera, scaleFactor=8):
    """Generate common plots of test failures.

    Parameters
    ----------
    runStats : `dict` [`str` `dict`]
        Statistics dictionary to retrieve failures from.
    camera : `lsst.afw.cameraGeom.Camera`
        Camera object to use for focal plane geometry.
    scaleFactor : `float`
        Scale factor to apply to plots.
    """
        numExposures = len(runStats.keys())
    # I probably should rethink the structure if I need to do this much work.
    failedTests = defaultdict(lambda: defaultdict(lambda: defaultdict(list))) # test -> detector -> amp -> list of exposure #
    for exposure, stats in runStats.items():
        failures = stats['FAILURES']
        for failure in failures:
            detector, amp, test = failure.split(" ")
            failedTests[test][detector][amp].append(exposure)

    # Focal plane plots
    numTests = len(failedTests.keys())
    fig, axes = plt.subplots(1, numTests, figsize=(numTests * scaleFactor, scaleFactor))
    if numTests == 1:
        axes = [axes]
    for axis, (testName, failureSet) in zip(axes, failedTests.items()):
        # full camera
        ccdList = list(set([camera[detName].getId() for detName in failureSet.keys()]))

        metricPerCcdDict = defaultdict(float)
        for detName, detFails in failureSet.items():
            detId = camera[detName].getId()
            for amp, expList in detFails.items():
                metricPerCcdDict[str(detId)] += len(expList)

        axis.set_aspect("equal")
        plotCameraOutline(axis, camera, list(metricPerCcdDict.keys()),
                          metricPerCcdDict=metricPerCcdDict, metricStr=testName, fig=fig, fontSize=10)

    fig, axes = plt.subplots(1, numTests, figsize=(numTests * scaleFactor, scaleFactor))
    if numTests == 1:
        axes = [axes]
    for axis, (testName, failureSet) in zip(axes, failedTests.items()):
        for detName, detFails in failureSet.items():
            detector = camera[detName]
            # Plot Nfailures / chip / amp
            ampNames = [amp.getName() for amp in detector]
            ampValues = [len(detFails.get(ampName, [])) for ampName in ampNames]
            bins = list(map(lambda x: x-0.5, range(1,len(ampNames)+1)))

            axis.bar(ampNames, ampValues, width=1.0)
            axis.axhline(numExposures)
            axis.set_title(f"{testName} {detName}")
            axis.set_xticklabels(ampNames, rotation=45, rotation_mode='anchor', ha='right')
            axis.grid()

