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
import matplotlib.pyplot as plt

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable

import lsst.geom
import lsst.afw.math as afwMath
from lsst.obs.lsst import LsstCam
from lsst.afw import cameraGeom

__all__ = [
    "det_to_index",
    "mergeStatDict",
    "get_amp_patches",
    "plot_amp_boundaries",
    "plot_detector",
    "nsigma_range",
    "get_median_nsigma_range",
    "get_median_nsigma_range",
    "plot_fp_statistic",
]

# All LSST cameras will have an R22. Therefore,
# ComCam and LATISS will both map to the center
# raft in the detailed plots.
det_to_index = {
    190: (2, 0),  # R00
    189: (0, 2),
    191: (2, 2),
    192: (1, 2),
    0: (0, 3),  # R01
    1: (0, 4),
    2: (0, 5),
    3: (1, 3),
    4: (1, 4),
    5: (1, 5),
    6: (2, 3),
    7: (2, 4),
    8: (2, 5),
    9: (0, 6),  # R02
    10: (0, 7),
    11: (0, 8),
    12: (1, 6),
    13: (1, 7),
    14: (1, 8),
    15: (2, 6),
    16: (2, 7),
    17: (2, 8),
    18: (0, 9),  # R03
    19: (0, 10),
    20: (0, 11),
    21: (1, 9),
    22: (1, 10),
    23: (1, 11),
    24: (2, 9),
    25: (2, 10),
    26: (2, 11),
    194: (0, 12),  # R04
    193: (2, 14),
    195: (2, 12),
    196: (2, 13),
    27: (3, 0),  # R10
    28: (3, 1),
    29: (3, 2),
    30: (4, 0),
    31: (4, 1),
    32: (4, 2),
    33: (5, 0),
    34: (5, 1),
    35: (5, 2),
    36: (3, 3),  # R11
    37: (3, 4),
    38: (3, 5),
    39: (4, 3),
    40: (4, 4),
    41: (4, 5),
    42: (5, 3),
    43: (5, 4),
    44: (5, 5),
    45: (3, 6),  # R12
    46: (3, 7),
    47: (3, 8),
    48: (4, 6),
    49: (4, 7),
    50: (4, 8),
    51: (5, 6),
    52: (5, 7),
    53: (5, 8),
    54: (3, 9),  # R13
    55: (3, 10),
    56: (3, 11),
    57: (4, 9),
    58: (4, 10),
    59: (4, 11),
    60: (5, 9),
    61: (5, 10),
    62: (5, 11),
    63: (3, 12),  # R14
    64: (3, 13),
    65: (3, 14),
    66: (4, 12),
    67: (4, 13),
    68: (4, 14),
    69: (5, 12),
    70: (5, 13),
    71: (5, 14),
    72: (6, 0),  # R20
    73: (6, 1),
    74: (6, 2),
    75: (7, 0),
    76: (7, 1),
    77: (7, 2),
    78: (8, 0),
    79: (8, 1),
    80: (8, 2),
    81: (6, 3),  # R21
    82: (6, 4),
    83: (6, 5),
    84: (7, 3),
    85: (7, 4),
    86: (7, 5),
    87: (8, 3),
    88: (8, 4),
    89: (8, 5),
    90: (6, 6),  # R22
    91: (6, 7),
    92: (6, 8),
    93: (7, 6),
    94: (7, 7),
    95: (7, 8),
    96: (8, 6),
    97: (8, 7),
    98: (8, 8),
    99: (6, 9),  # R23
    100: (6, 10),
    101: (6, 11),
    102: (7, 9),
    103: (7, 10),
    104: (7, 11),
    105: (8, 9),
    106: (8, 10),
    107: (8, 11),
    108: (6, 12),  # R24
    109: (6, 13),
    110: (6, 14),
    111: (7, 12),
    112: (7, 13),
    113: (7, 14),
    114: (8, 12),
    115: (8, 13),
    116: (8, 14),
    117: (9, 0),  # R30
    118: (9, 1),
    119: (9, 2),
    120: (10, 0),
    121: (10, 1),
    122: (10, 2),
    123: (11, 0),
    124: (11, 1),
    125: (11, 2),
    126: (9, 3),  # R31
    127: (9, 4),
    128: (9, 5),
    129: (10, 3),
    130: (10, 4),
    131: (10, 5),
    132: (11, 3),
    133: (11, 4),
    134: (11, 5),
    135: (9, 6),  # R32
    136: (9, 7),
    137: (9, 8),
    138: (10, 6),
    139: (10, 7),
    140: (10, 8),
    141: (11, 6),
    142: (11, 7),
    143: (11, 8),
    144: (9, 9),  # R33
    145: (9, 10),
    146: (9, 11),
    147: (10, 9),
    148: (10, 10),
    149: (10, 11),
    150: (11, 9),
    151: (11, 10),
    152: (11, 11),
    153: (9, 12),  # R34
    154: (9, 13),
    155: (9, 14),
    156: (10, 12),
    157: (10, 13),
    158: (10, 14),
    159: (11, 12),
    160: (11, 13),
    161: (11, 14),
    197: (12, 0),  # R40
    198: (14, 2),
    199: (12, 2),
    200: (12, 1),
    162: (12, 3),  # R41
    163: (12, 4),
    164: (12, 5),
    165: (13, 3),
    166: (13, 4),
    167: (13, 5),
    168: (14, 3),
    169: (14, 4),
    170: (14, 5),
    171: (12, 6),  # R42
    172: (12, 7),
    173: (12, 8),
    174: (13, 6),
    175: (13, 7),
    176: (13, 8),
    177: (14, 6),
    178: (14, 7),
    179: (14, 8),
    180: (12, 9),  # R43
    181: (12, 10),
    182: (12, 11),
    183: (13, 9),
    184: (13, 10),
    185: (13, 11),
    186: (14, 9),
    187: (14, 10),
    188: (14, 11),
    201: (14, 12),  # R44
    202: (12, 14),
    203: (12, 12),
    204: (13, 12),
}


def mergeStatDict(statDictA, statDictB):
    """Merge one dictionary of statistics with another.

    Parameters
    ----------
    statDictA, statDictB: `dict` [`str`, `dict`]
        Input dictionaries of statistics.

    Returns
    -------
    output: `dict` [`str`, `dict`]
        A dictionary containing the union of the internal
        dictionaries of ``statDictA`` and ``statDictB``, indexed
        by the amplifier names.

    Raises
    ------
    RuntimeError:
        Raised if the merge would overwrite a key/value pair, or if
        the amplifier names do not match between the two inputs.
    """
    if set(statDictA.keys()) != set(statDictB.keys()):
        raise RuntimeError(f"Mismatch in amplifier names: {statDictA.keys()} {statDictB.keys()}")

    for amp in statDictA.keys():
        if not statDictA[amp].keys().isdisjoint(statDictB[amp].keys()):
            raise RuntimeError(f"Duplicate keys passed: {statDictA[amp].keys()} {statDictB[amp].keys()}")
        statDictA[amp].update(statDictB[amp])
    return statDictA


# Plotting
def get_amp_patches(det, amps=None):
    """
    Return a list of Rectangle patches in focalplane coordinates
    corresponding to the amplifier segments in a detector object.

    Parameters
    ----------
    det: `lsst.afw.cameraGeom.detector.detector.Detector`
        Detector object.
    amps: container-type object [None]
        Python container that can be queried like `'C01 in amps'`
        to see if a particular channel is included for plotting.
        If None, then use all channels in det.

    Returns
    -------
    list of matplotlib.patches.Rectangle objects
    """
    transform = det.getTransform(cameraGeom.PIXELS, cameraGeom.FOCAL_PLANE)
    bbox = list(det)[0].getBBox()
    dy, dx = bbox.getHeight(), bbox.getWidth()
    patches = []
    for amp in det:
        if amps is not None and amp.getName() not in amps:
            continue
        j, i = tuple(int(_) for _ in amp.getName()[1:])
        y, x = j*dy, i*dx
        x0, y0 = transform.applyForward(lsst.geom.Point2D(x, y))
        x1, y1 = transform.applyForward(lsst.geom.Point2D(x + dx, y + dy))
        patches.append(Rectangle((x0, y0), x1 - x0, y1 - y0))
    return patches


def plot_amp_boundaries(ax, camera=None, edgecolor='gray', facecolor='white'):
    """
    Plot the amplifier boundaries for all of the detectors in a camera.

    Parameters
    ----------
    ax: `matplotlib.Axes`
        Axes object used to render the patch collection containing
        the amplifier boundary Rectangles.
    camera: `lsst.afw.cameraGeom.camera.camera.Camera` [None]
        Camera object containing the detector info. If None, use
        `lsst.obs.lsst.LsstCam.getCamera()`
    edgecolor: str or tuple of RGBA values ["blue"]
        Color used to draw outline of amplifiers.
    facecolor: str or tuple of RGBA values ["white"]
        Color used to render the Rectangle corresponding to the
        amplifier region.

    Returns
    -------
    None
    """
    if camera is None:
        camera = LsstCam.getCamera()
    patches = []
    for det in camera:
        patches.extend(get_amp_patches(det))
    pc = PatchCollection(patches, edgecolor=edgecolor, facecolor=facecolor)
    ax.add_collection(pc)


def plot_detector(ax, det, amp_values, cm=plt.cm.hot, z_range=None, use_log10=False):
    """
    Plot the amplifiers in a detector, rendering each amplier region with
    a color corresponding to its assigned value.

    Parameters
    ----------
    ax: `matplotlib.Axes`
        Axes object used to render the patch collection containing
        the amplifier boundary Rectangles.
    det: `lsst.afw.cameraGeom.detector.detector.Detector`
        Detector object.
    amp_values: dict of floats
        Dictionary of amplifier values to render, keyed by channel ID,
        e.g., 'C00', 'C01', etc.
    cm: `matplotlib.colors.Colormap`
        Colormap used to render amplifier values.
    z_range: 2-tuple of floats [None]
        Minimum and maximum values into which to map the unit interval
        for the color map.  Values are mapped using
        max(0, min(1, (amp_value - z_range[0])/(z_range[1] - z_range[0])))
        If None, then use
        z_range = (min(amp_values.values()), max(amp_values.values()))
    use_log10: bool [False]
        If True, then use log10(amp_value) for positive amp_value.  For
        non-positive amp_value, don't render the amp color.

    Returns
    -------
    None
    """
    if z_range is None:
        zvals = amp_values.values()
        z_range = min(zvals), max(zvals)

    def mapped_value(amp_value):
        return max(0, min(1., ((amp_value - z_range[0])
                               / (z_range[1] - z_range[0]))))

    my_facecolors = []
    for amp in det:
        if amp.getName() not in amp_values:
            continue
        if use_log10:
            if amp_values[amp.getName()] <= 0:
                my_facecolors.append(None)
            else:
                my_facecolors.append(
                    cm(mapped_value(np.log10(amp_values[amp.getName()]))))
        else:
            my_facecolors.append(cm(mapped_value(amp_values[amp.getName()])))

    my_patches = get_amp_patches(det, amp_values)
    facecolors, patches = [], []
    for facecolor, patch in zip(my_facecolors, my_patches):
        if facecolor is not None:
            facecolors.append(facecolor)
            patches.append(patch)
    assert len(facecolors) == len(patches)
    pc = PatchCollection(patches, facecolors=facecolors)
    ax.add_collection(pc)


def nsigma_range(data, nsigma=3):
    stats = afwMath.makeStatistics(
        data,
        afwMath.MEDIAN | afwMath.STDEVCLIP | afwMath.STDEV,
    )
    median = stats.getValue(afwMath.MEDIAN)
    stdev = stats.getValue(afwMath.STDEVCLIP)
    if not np.isfinite(stdev):
        return None
    return (median - nsigma*stdev, median + nsigma*stdev)


def get_median_nsigma_range(amp_data, nsigma=4, use_log10=False):
    """
    Use a clipped stdev to compute the range of amp_data appropriate
    for plotting as a full focal plane mosaic or as a histogram of
    per-amp values.  The returned range will be
    (median - nsigma*stdev_clip, median + nsigma*stdev_clip).

    Parameters
    ----------
    amp_data: dict of dict of floats
        Dictionary of dictionary of amplifier values to render,
        keyed by detector name, e.g., 'R01_S00' and then by channel ID,
        e.g., 'C00'.
    nsigma: float [4]
        Number of sigma to use for computing the plotting range.
    use_log10: bool [False]
        If True, then use the log10 of the amp_data values, excluding
        non-positive values.

    Returns
    -------
    (float, float)
    """
    amp_values = []
    for _ in amp_data.values():
        amp_values.extend(_.values())
    amp_values = np.array(amp_values, dtype=np.float64)
    if use_log10:
        amp_values = np.array([np.log10(_) for _ in amp_values if _ > 0])
    return nsigma_range(amp_values, nsigma=nsigma)


def plot_fp_statistic(ax, amp_data, camera=None, cm=plt.cm.hot,
                      x_range=None, y_range=None,
                      z_range=None, use_log10=False, scale_factor='1',
                      title='', nsigma=4):
    """
    Make a "heat map" plot of the focalplane using per-amplifier values.

    Parameters
    ----------
    ax: `matplotlib.Axes`
        Axes object used to render the patch collection containing
        the amplifier boundary Rectangles.
    amp_data: dict of dict of floats
        Dictionary of dictionary of amplifier values to render,
        keyed by detector name, e.g., 'R01_S00' and then by channel ID,
        e.g., 'C00'.
    camera: `lsst.afw.cameraGeom.camera.camera.Camera` [None]
        Camera object containing the detector info. If None, use
        `lsst.obs.lsst.LsstCam.getCamera()`
    cm: `matplotlib.colors.Colormap`
        Colormap used to render amplifier values.
    x_range: tuple [None]
        Focalplane plotting region in x-direction in units of mm.
        If None, then the region will be inferred from the camera.
    y_range: tuple [None]
        Focalplane plotting region in y-direction in units of mm.
        If None, then the region will be inferred from the camera.
    z_range: 2-tuple of floats or `str` [None]
        Minimum and maximum values into which to map the unit interval
        for the color map.  Values are mapped using
        max(0, min(1, (amp_value - z_range[0])/(z_range[1] - z_range[0])))
        If None, then use
        z_range = (min(amp_data.values(), max(amp_data.values())
        If "clipped_autoscale", then use
        z_range = get_median_nsigma_range(amp_data, nsigma=nsigma,
                                            use_log10=use_log10)
    use_log10: bool [False]
        If True, then use log10(amp_value) for positive amp_value.  For
        non-positive amp_value, don't render the amp color.
    scale_factor: str ['1']
        Scale factor to apply to the colorbar mapping.  This value
        will be cast as a float when applied to the tick label values.  It
        is passed as a string so that formatting in the colorbar tick
        labels can be controlled directly by the client code.
        This is not used if use_log10 == True.
    title: str ['']
        Title to apply to the plot.
    nsigma: float [4]
        Number of sigma to apply to clipped stdev for an "autoscaled"
        z_range, i.e., [median - nsigm*stdev_clip, median + nsigm*stdev_clip].

    Returns
    -------
    None
    """
    xy_ranges = {'LSSTCam': (-325, 325),
                 'LSST-TS8': (-70, 70)}
    if camera is None:
        camera = LsstCam.getCamera()
    fp_radius = camera.computeMaxFocalPlaneRadius()
    if x_range is None:
        x_range = xy_ranges.get(camera.getName(), (-fp_radius, fp_radius))
    if y_range is None:
        y_range = xy_ranges.get(camera.getName(), (-fp_radius, fp_radius))
    plot_amp_boundaries(ax, camera=camera)
    if z_range is None:
        z_range_values = []
        for _ in amp_data.values():
            z_range_values.extend(_.values())
        if use_log10:
            z_range_values = [np.log10(_) for _ in z_range_values if _ > 0]
        z_range = min(z_range_values), max(z_range_values)
    elif z_range == 'clipped_autoscale':
        z_range = get_median_nsigma_range(amp_data, nsigma=nsigma,
                                          use_log10=use_log10)
    for det_name, amp_values in amp_data.items():
        plot_detector(ax, camera[det_name], amp_values, cm=cm, z_range=z_range,
                      use_log10=use_log10)
    ax.set_xlim(*x_range)
    ax.set_ylim(*y_range)
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    norm = plt.Normalize(vmin=z_range[0], vmax=z_range[1])
    sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
    sm.set_array([])
    ax.set_title(title)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    colorbar = plt.colorbar(sm, cax=cax)
    if use_log10:
        ticks = sorted(list(set([int(_) for _ in
                                np.log10(np.logspace(0, z_range[1]))])))
        ticklabels = [10**_ for _ in ticks]
        colorbar.set_ticks(ticks)
        colorbar.set_ticklabels(ticklabels)
    elif scale_factor != '1':
        ticks = colorbar.get_ticks()
        colorbar.set_ticks(ticks)
        ticklabels = [_/float(scale_factor) for _ in ticks]
        ticklabels[-1] = '{} x {}'.format(ticklabels[-1], scale_factor)
        colorbar.set_ticklabels(ticklabels)
    ax.set_aspect('equal')
    return colorbar
