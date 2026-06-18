# This file is part of cp_verify.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This software is dual licensed under the GNU General Public License and also
# under a 3-clause BSD license. Recipients may choose which of these licenses
# to use; please see the files gpl-3.0.txt and/or bsd_license.txt,
# respectively.  If you choose the GPL option then the following text applies
# (but note that there is still no warranty even if you opt for BSD instead):
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
import os
import re
import argparse
import logging

from collections import defaultdict

from lsst.daf.butler import Butler
from lsst.obs.lsst import LsstCam

from .utils import det_to_index


class DetailedReporterLSSTCam():
    """A class to generate calibration verification reports.

    Parameters
    ----------
    sharedPath : `str`
        The location the report will be written to.
    instrument : `str`
        The instrument associated with the data.
    doOverwrite : `bool`
        Should the sharedPath be overwritten?
    """

    def __init__(self, sharedPath, instrument, doOverwrite):
        super().__init__()
        # Set source and destination information.
        self.sharedPath = sharedPath
        self.doOverwrite = doOverwrite

        # Check that the output directory
        # is organized in the proper format
        if not os.path.isdir(self.sharedPath):
            raise RuntimeError(f"The shared directory {self.sharedPath} does not exist.")
        if not os.path.isdir(self.sharedPath + "/images"):
            raise RuntimeError(f"The shared directory {self.sharedPath} does not follow the expected "
                               "organizational structure, please refer to DMTN-222 for more information.")

        # Logging
        self.log = logging.getLogger(__name__)

        # Instantiate a butler for our repository.
        self.butler = Butler("/repo/main")
        self.registry = self.butler.registry

        # Get the camera associated with these calibrations
        if instrument != "LSSTCam":
            raise RuntimeError("The only supported camera is currently LSSTCam.")
        self.camera = LsstCam.getCamera()

    def writeReport(self, renderedHtml, fileName):
        output_path = self.sharedPath + "/" + fileName

        # Do nothing if the file already exists and doOverwrite==True
        if os.path.isfile(output_path) and not self.doOverwrite:
            self.log.warning(f"The file {output_path} already exists and doOverwrite=False, skipping.")
            return

        with open(output_path, "w") as f:
            f.write(renderedHtml)

    def getAllImages(self):
        summaryImages = self.getSummaryImages()
        detailedImages = self.getDetailedImages()

        return summaryImages, detailedImages

    def getSummaryImages(self):
        mosaicImages = []
        histImages = []

        # Get focal plane images
        imagePath = os.path.join(self.sharedPath, "images")
        if os.path.exists(imagePath):
            for f in os.listdir(imagePath):
                if f.endswith('hist_mosaic.png') and os.path.isfile(os.path.join(imagePath, f)):
                    histImages.append(f)
                elif f.endswith('_mosaic.png') and os.path.isfile(os.path.join(imagePath, f)):
                    mosaicImages.append(f)
                else:
                    continue

        if len(mosaicImages) != len(histImages):
            raise RuntimeError("The number of summary mosaic images does not match "
                               "the number of histogram images.")

        mosaicImages = sorted(mosaicImages)
        histImages = sorted(histImages)

        summaryImages = list(zip(mosaicImages, histImages))

        return summaryImages

    def getDetailedImages(self):
        tree = defaultdict(list)
        # Get tree data
        fullNames = []
        for root, dirs, files in os.walk(self.sharedPath):
            for file in files:
                if file.endswith('.png') and not file.endswith('_fp.png'):
                    try:
                        pattern = re.compile(r"(step_\d{2}[a-zA-Z]?)-([a-zA-Z0-9_-]+)_(\d{2,3})\.png")
                        match = pattern.match(file)
                        if match is None:
                            continue
                        step = match.group(1)
                        name = match.group(2)
                        fullName = f"{step}-{name}"
                        fullNames.append(fullName)
                        detId = int(match.group(3))
                        relPath = os.path.join(".", os.path.relpath(root, self.sharedPath), file)
                        tree[fullName].append((detId, relPath))
                    except (ValueError, IndexError):
                        continue

        # Sort tree data
        for fullName in tree:
            tree[fullName].sort(key=lambda x: x[0])
        tree = dict(sorted(tree.items()))

        # Process tree data into grids
        detailedImages = {}
        for fullName, images in tree.items():
            # Create 15x15 grid (initialized with None)
            grid = [[None for _ in range(15)] for _ in range(15)]
            for detId, path in images:
                if detId in det_to_index:
                    r, c = det_to_index[detId]
                    grid[14 - r][c] = path
            detailedImages[fullName] = grid

        return detailedImages

    def makeSummary(self, summaryImages):
        # Get the template for the summary report
        summaryReportTemplateFile = os.path.join(
            os.environ['CP_VERIFY_DIR'],
            "python/lsst/cp/verify/templates/summary.html",
        )
        with open(summaryReportTemplateFile, 'r') as file:
            summaryReportTemplate = file.read()

        addition = ""
        for i, (mosaicPath, histPath) in enumerate(summaryImages):
            mosaicPath = "./images/" + mosaicPath
            histPath = "./images/" + histPath
            addition += """
                <a href="{0}" target="_blank">
                    <img src="{0}" alt="Mosaic Image {2}">
                </a>
                <a href="{1}" target="_blank">
                    <img src="{1}" alt="Hist Image {2}">
                </a>
            """.format(mosaicPath, histPath, i)

        renderedHtml = summaryReportTemplate.replace("summary_images", addition)

        self.writeReport(renderedHtml, "summary.html")

    def makeDetailed(self, fullName, tree):
        """A class to generate the detailed verification reports.

        Parameters
        ----------
        fullName : `str`
            The name of the file to output.

        tree :  `dict` [`list`, `tuple`]
            Dictionary keyed by "{step}-{name}", that contains
            a list of tuples (detectorId, relative/path/to/image).

        """
        # Get the template for the detailed report
        detailedReportTemplateFile = os.path.join(
            os.environ['CP_VERIFY_DIR'],
            "python/lsst/cp/verify/templates/detailed.html",
        )
        with open(detailedReportTemplateFile, 'r') as file:
            detailedReportTemplate = file.read()

        pattern = re.compile(r"^step_\d{2}[a-zA-Z]?-(.+)$")
        match = pattern.match(fullName)
        name = match.group(1)
        name = name.replace("_", " ").upper()
        name = name.replace("AND", "VS").upper()

        addition = ""
        cell = """
        <div class="grid-cell">
                cell_image
        </div>
        """
        for i, paths in enumerate(tree):
            for path in paths:
                if path is not None:
                    img = """
                    <a href="{0}" target="_blank">
                        <img src="{0}"
                             alt="{1} image"
                             title="Click to view full size">
                    </a>
                    """.format(path, "test")
                    addition += cell.replace("cell_image", img)
                else:
                    addition += cell.replace("cell_image", "")
                    continue

        renderedHtml = detailedReportTemplate.replace("grid_placeholder", addition).replace("name", name)

        self.writeReport(renderedHtml, f"{fullName}.html")

    def run(self):
        "Get all the images and generate the HTML construction routine."
        # What images are there?
        summaryImages, detailedImages = self.getAllImages()

        # Make summary plots
        if len(summaryImages) > 0:
            self.makeSummary(summaryImages)

        # Make detailed plots
        if len(detailedImages) > 0:
            for fullName in detailedImages.keys():
                tree = detailedImages[fullName]
                self.makeDetailed(fullName, tree)


def main():
    r"""
    Construct a webpage (html) report of a full focal plane given
    a set of images for each detector.
    The reporter only requires that the images are organized in a
    particular format (the images directory can be anywhere, but
    preferably in ones public_html directory on S3D)F:

    images/
        |_____R00/
        |          |______S00/
        |          |______S01/
        |               ...
        |_____R01/ ...
        |_____ ...

    Also, the images in the directory need to follow this naming
    convention:

    ``step_\d{2}[a-zA-Z]?)-([a-zA-Z0-9_-]+)_(\d{2,3})\.png``
    ==
    step_(## + letter)-(descriptive name)_(detector ###).png

    Or, you can adjust line 132 on your local version of
    python/lsst/cp/verify/detailedReportsLSSTCam.py to fit your
    own format.

    The report code will generate a .html file for each unique
    step_(## + letter)-(descriptive name) name.
    """
    parser = argparse.ArgumentParser(
        description="Construct a detailed report. Usage: -p (absolute path "
                    "to output directory) --instrument (e.g. LSSTCam) "
                    "--do-overwrite (overwrite the path where the files are "
                    "written.)")
    parser.add_argument("-p", "--shared-path", dest="sharedPath", default="",
                        help="Absolute path with input 'images' diectory, which will also "
                             "contain the report pages.")
    parser.add_argument("--instrument", dest="instrument", default="",
                        help="The instrument to use (e.g. LSSTCam)")
    parser.add_argument("--do-overwrite", dest="doOverwrite", action="store_true", default=False,
                        help="Allow existing files to be overwritten?")
    args = parser.parse_args()

    reporter = DetailedReporterLSSTCam(
        sharedPath=args.sharedPath,
        instrument=args.instrument,
        doOverwrite=args.doOverwrite,
    )
    reporter.run()
