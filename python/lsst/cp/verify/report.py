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
import argparse
import logging
import matplotlib.colors as colors
import numpy as np
import os
import re
import yaml

from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from lsst.daf.butler import Butler
from lsst.daf.butler.datastores.file_datastore.retrieve_artifacts import (
    determine_destination_for_retrieved_artifact,
)
from lsst.resources import ResourcePath
from lsst.utils import getPackageDir
from lsst.utils.plotting import make_figure


class CpvReporter():
    """A class to generate calibration verification reports.

    Parameters
    ----------
    repo : `str`
        The location of the butler repository to retrieve results
        from.
    instrument : `str`
        The instrument associated with this data.
    output_path : `str`
        The location the report will be written to.
    collections : `list` [`str`]
        A list of collections to search.
    **kwargs :
        Other keyword parameters.  Currently parsed values:

        - ``do_copy``: Should files be copied from butler (`bool`).
        - ``do_overwrite``: Should pre-existing files be overwritten (`bool`).
    """

    def __init__(self, repo, output_path, collections=[], **kwargs):
        super().__init__()
        # Set source and destination information.
        self.repo = repo
        self.output_path = output_path
        self.src_dir = os.path.join(self.output_path, "src")
        self.collections = collections

        # We need a log
        self.log = logging.getLogger(__name__) if "log" not in kwargs else kwargs["log"]
        # List of the datasets we need to process.  The datasets are
        # dictionaries of parameters, including the dataset_type_name,
        # the source collection, the butler dataIds, file locations,
        # etc.
        self.datasets = []

        # A dictionary of output html files that will be created.
        # These are dictionaries with the key being the filename, and
        # the values being arrays of strings that will be written
        # one-per-line in the output html page.
        self.out_files = {
            'index.html': self._init_page(),
        }

        # Instantiate a butler for our repository.
        self.butler = Butler(repo)

        # Configure behavior from kwargs:
        self.do_copy = kwargs['do_copy'] if 'do_copy' in kwargs else True
        self.do_overwrite = kwargs['do_overwrite'] if 'do_overwrite' in kwargs else False

        # Read dataset configuration yaml:
        self.dataset_map = self._read_dataset_map()

    def _read_dataset_map(self):
        """Read dataset information from source yaml."""
        filename = os.path.join(getPackageDir("cp_verify"),
                                "python", "lsst", "cp", "verify", "configs", "report.yaml")
        with open(filename) as in_file:
            return yaml.safe_load(in_file)

    def run(self):
        """Generate the report"""
        # Generate directories that will hold the output products.
        os.makedirs(self.src_dir, exist_ok=True)

        # Copy datasets:  This populates self.datasets.
        self.copy_datasets()

        # Parse datasets:  This determines what report pages should be made.
        self.parse_datasets()

        # Write pages:  This saves the information to disk.
        self.write_pages()

    def copy_datasets(self):
        """Copy datasets, similar to butler retrieve-artifacts.

        See Also
        --------
        self.fits_to_png
        self.parq_to_html
        self.py_to_html
        """
        # Based on daf_butler fileDatastore.py#L1978
        for stage, dataset_types in self.dataset_map['stages'].items():
            # Iterate over each stage, getting all of that stage's datasets.
            # Convert dict dataset_types to a list of strings:
            if not dataset_types:
                continue
            dataset_type_names = list(dataset_types.keys())

            # Find all the references to these datasets in our collections
            refs = self.butler.registry.queryDatasets(datasetType=dataset_type_names,
                                                      collections=self.collections,
                                                      where=None, findFirst=True)

            for ref in refs:
                # Iterate over each reference
                locations = self.butler._datastore._get_dataset_locations_info(ref)
                self.log.warn(f"Found {stage} {ref}")

                # This is our dataset information:
                dataset = {'stage': stage,
                           'ref': ref,
                           'type': ref.datasetType.name,
                           'storage_class': ref.datasetType.storageClass.name,
                           'dataId': ref.dataId,
                           'instrument': ref.dataId.get('instrument', None),
                           'exposure': ref.dataId.get('exposure', None),
                           'detector': ref.dataId.get('detector', None),
                           'physical_filter': ref.dataId.get('physical_filter', None),
                           'collection': ref.run,
                           }

                # Do the actual copy.  This may convert the butler
                # product to something more web-accessible.
                for location, _ in locations:
                    source_uri = location.uri
                    target_uri = determine_destination_for_retrieved_artifact(
                        ResourcePath(self.src_dir),
                        location.pathInStore,
                        False
                    )
                    dataset['uri'] = target_uri
                    willCopy = self.do_copy

                    if ((target_uri.getExtension().lower() == ".fits"
                         and dataset['storage_class'] in ('ImageF', 'ExposureF', 'ImageD', 'ExposureD'))):
                        # Convert this to a PNG image.
                        png_uri = target_uri.updatedExtension("png")
                        dataset['uri'] = png_uri
                        if os.path.exists(png_uri.path):
                            if willCopy and not self.do_overwrite:
                                willCopy = False

                        if willCopy:
                            data = self.butler.get(ref)
                            self.fits_to_png(data, png_uri, f"{ref.datasetType.name} {ref.dataId}")
                    elif (target_uri.getExtension().lower() == ".parq"
                          and dataset['storage_class'] in ('ArrowAstropy', )):
                        # Convert to html table.
                        html_uri = target_uri.updatedExtension("html")
                        dataset['uri'] = html_uri
                        if os.path.exists(html_uri.path):
                            if willCopy and not self.do_overwrite:
                                willCopy = False

                        if willCopy:
                            data = self.butler.get(ref)
                            self.parq_to_html(data, html_uri)
                    elif (target_uri.getExtension().lower() == ".py"
                          and dataset['storage_class'] in ('Config', )):
                        # Convert configs to html.
                        html_uri = target_uri.updatedExtension("html")
                        dataset['uri'] = html_uri

                        if willCopy:
                            data = self.butler.get(ref)
                            self.py_to_html(data, html_uri)
                    else:
                        if willCopy:
                            # Otherwise just copy directly.
                            target_uri.transfer_from(source_uri, transfer="copy",
                                                     overwrite=self.do_overwrite)

                # Add this dataset to our list.
                self.datasets.append(dataset)

    def parse_datasets(self):
        """Analyze datasets, and add information to report pages.

        See Also
        --------
        self._init_page
        self.include
        self.navigation_init
        self.navigation_append
        self.title
        self.link
        """
        populated_stages_list = []

        # Make navigation page.
        self.out_files["navigation.html"] = self._init_page("Navigation")
        self.navigation_init()
        self.include("navigation.html", "index.html")

        # Make manifest page:
        self.out_files["manifest.html"] = self._init_page("Manifest")

        tag = self.link(self.out_files["index.html"], "Parent", "./..")
        self.navigation_append("index.html", "Parent", tag, newColumn=True)
        self.navigation_append("manifest.html", "Manifest", "", newColumn=False)

        for stage in self.dataset_map['stages'].keys():
            # Filter datasets to only those for this stage
            dss = [x for x in self.datasets if x['stage'] == stage]
            if len(dss) > 0:
                populated_stages_list.append(stage)
            else:
                continue
            self.log.warn(f"Parsing data for {stage}")
            self.title(self.out_files["index.html"], stage)

            self.out_files[f"{stage}_exp.html"] = self._init_page()
            self.include("navigation.html", f"{stage}_exp.html")
            self.out_files[f"{stage}_det.html"] = self._init_page()
            self.include("navigation.html", f"{stage}_det.html")

            # Get dimension information we'll want to use to put
            # things on appropriate pages.
            exps = [x for x in dss if x['exposure'] is not None]
            dets = [x for x in dss if x['detector'] is not None]

            for ds in sorted(exps, key=lambda x: x['exposure']):
                self.block(ds, self.out_files[f"{stage}_exp.html"])
                if ('.png' in ds['uri'].path):
                    self.image_handler(ds['uri'].path, ds['instrument'], self.out_files[f"{stage}_exp.html"])
            tag = self.link(self.out_files["index.html"], f"{stage} exposures", f"./{stage}_exp.html")
            self.navigation_append("index.html", f"{stage} exposures", tag)

            for ds in sorted(dets, key=lambda x: x['detector']):
                self.block(ds, self.out_files[f"{stage}_det.html"])
                if ('.png' in ds['uri'].path):
                    self.image_handler(ds['uri'].path, ds['instrument'], self.out_files[f"{stage}_det.html"])
            tag = self.link(self.out_files["index.html"], f"{stage} detectors", f"./{stage}_det.html")
            self.navigation_append("index.html", f"{stage} detectors", tag)

            for ds in dss:
                if ds['storage_class'] == 'ArrowAstropy':
                    # Is a results catalog
                    self.block(ds, self.out_files["index.html"], ds['uri'])
                    self.include(ds['uri'].path, "index.html")
                elif ds['exposure'] is None and ds['detector'] is None:
                    self.block(ds, self.out_files["index.html"])
                    if ".png" in ds['uri'].path:
                        self.image_handler(ds['uri'].path, ds['instrument'], self.out_files["index.html"])
                    else:
                        self.include(ds['uri'].path, "index.html")

        # Record all datasets here, and link to them directly:
        page = self.out_files["manifest.html"]
        page.extend(["<table width='100%'>",
                     "<tr>",
                     "<th>stage</th>", "<th>ref</th>", "<th>type</th>", "<th>storage_class</th>",
                     "<th>dataId</th>", "<th>collection</th>",
                     "</tr>"])
        for ds in self.datasets:
            relative_file = self.relative_file(ds['uri'].path)
            page.extend(["<tr>",
                         f"<td>{ds['stage']}</td>", f"<td>{ds['ref']}</td>",
                         "<td>", f'<a href="{relative_file}">', f"{ds['type']}", "</a>",
                         f"<td>{ds['storage_class']}</td>",
                         f"<td>{ds['dataId']}</td>", f"<td>{ds['collection']}</td>",
                         "</tr>"])

    # File conversion utilities:
    @staticmethod
    def fits_to_png(data, target_uri, title):
        """Convert FITS to PNG, using the same scaling as on RubinTV.

        Parameters
        ----------
        data : `lsst.afw.Image`, `lsst.afw.MaskedImage` or `lsst.afw.Exposure`
            The fits image data to convert.
        target_uri : `str`
            Path to the PNG file to write.
        title : `str`
            Title to add to the figure.

        See Also
        --------
        rubintv_production / mosaicing.py
        """
        # Get array from either an Image or an Exposure:
        try:
            array = data.image.array
        except AttributeError:
            array = data.array

        fig = make_figure()
        ax = fig.gca()
        ax.clear()
        cmap = cm.gray
        # This was using summit_utils.getQuantiles, but that
        # was blowing out the scaling more than I wanted.
        q25, q50, q75 = np.nanpercentile(array, [25, 50, 75])
        scale = 3.0 * 0.74 * (q75 - q25)
        quantiles = np.arange(q50 - scale, q50 + scale, 2.0 * scale / cmap.N)
        norm = colors.BoundaryNorm(quantiles, cmap.N)

        im = ax.imshow(array, norm=norm, interpolation='None', cmap=cmap, origin='lower')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.suptitle(title)
        fig.colorbar(im, cax=cax)
        fig.tight_layout()
        fig.savefig(target_uri.path)

    @staticmethod
    def parq_to_html(data, target_uri):
        """Convert catalogs to html tables.

        Parameters
        ----------
        data : `astropy.Table`
            The table to convert to html.
        target_uri : `str`
            Path to the HTML file to write.
        """
        format_dict = {}
        # Remove vectors, as they will be too big.
        # Build up format dictionary.
        columns_to_remove = []
        for index, name in enumerate(data.dtype.names):
            if len(data.dtype[index].shape) != 0:
                columns_to_remove.append(name)
                continue
            if data.dtype[index].kind in ('f', 'c'):
                # is float
                format_dict[name] = "%.4g"
                if name == 'mjd':
                    # Let's let these be long
                    format_dict[name] = "12.10f"
            elif data.dtype[index].kind in ('i', 'u'):
                # is int
                format_dict[name] = "%d"
            else:
                format_dict[name] = "%s"

        # Actually remove the vectors:
        data.remove_columns(columns_to_remove)

        # Write full table:
        data.write(target_uri.path,
                   format='ascii.html',
                   overwrite=True,
                   formats=format_dict)

        # TODO: Write subset tables, filtered by exposure and by detector.
        # files_created[f"{filename_base}_exp{exp}.html"] = {'exposure': exp}
        # files_created[f"{filename_base}_det{det}.html"] = {'detector': det}

    def py_to_html(self, data, target_uri):
        """Convert python configs to html.

        Parameters
        ----------
        data : `Config`
            The python file contents to write.
        target_uri : `str`
            Path to the HTML file to write.
        """
        with open(target_uri.path, 'w') as ff:
            for line in self._init_page():
                print(line, file=ff)
            for line in data.saveToString().split("\n"):
                if len(line) > 0 and line[0] == '#':
                    print("<pre class='pre-comment'>", line, "</pre>", file=ff)
                else:
                    print("<pre>", line, "</pre>", file=ff)
            for line in self._close_page():
                print(line, file=ff)

    def write_pages(self):
        """Write all page arrays to files."""
        for ff, contents in self.out_files.items():
            reportOut = os.path.join(self.output_path, ff)
            contents.extend(self._close_page())
            with open(reportOut, 'w') as ff:
                for line in contents:
                    print(line, file=ff)

    # Page construction methods.  These all return arrays of strings
    # that can be appended to until we write to disk.
    @staticmethod
    def _init_page(title=""):
        """Write boilerplate html for the start of a page."""
        return [
            "<html>", "<head>", "<style> ",
            " * { margin: 0; padding: 0;}",
            " .imgbox { display: grid; height: 100%; }",
            " .pre-comment { color: red; }",
            " td { padding: 0 15px; }",
            " .ctable { display: table; }",
            " .ctable tr { display: table-cell; }",
            " .ctable tr td { display: block; }",
            "</style>",
            "<title>", title, "</title>",
            '<base target="_parent">',
            "</head>",
            "<body>",
        ]

    @staticmethod
    def _close_page():
        """Write boilerplate html for the end of a page."""
        return ["</body>",
                "</html>"]

    @staticmethod
    def relative_file(filename):
        """Convert absolute paths to the relative format needed for the html.
        """
        return re.sub(r"^.*/src/", "./src/", filename)

    def svg_atools_latiss(self, relative_file, page):
        """This adds an image to the html, with an SVG overlay for LATISS.

        Parameters
        ----------
        relative_file : `str`
            Relative location of the image to be added.
        page : `list` [`str`]
            The page to append to.

        Notes
        -----
        This method assumes the image is made by analysis_tools, and
        displays a focal plane plot for LATISS.  Under these
        assumptions, the boxes defined in SVG map to the amplifier
        segments of the single detector.
        """
        page.extend(['<svg version="1.1" xmlns="http://www.w3.org/2000/svg" '
                     'xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 3600 1800">',
                     f'<image width="3600" height="1800" xlink:href="{relative_file}"></image>',
                     '<a xlink:href="./c10.html"> <rect x="317" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c11.html"> <rect x="490" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c12.html"> <rect x="670" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c13.html"> <rect x="845" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c14.html"> <rect x="1020" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c15.html"> <rect x="1197" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c16.html"> <rect x="1375" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c17.html"> <rect x="1550" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c07.html"> <rect x="1550" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c06.html"> <rect x="1375" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c05.html"> <rect x="1197" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c04.html"> <rect x="1020" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c03.html"> <rect x="845" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c02.html"> <rect x="670" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c01.html"> <rect x="490" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c00.html"> <rect x="317" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '</svg>'])

    def svg_atools_comcam(self, relative_file, page):
        """This adds an image to the html, with an SVG overlay for ComCam.

        Parameters
        ----------
        relative_file : `str`
            Relative location of the image to be added.
        page : `list` [`str`]
            The page to append to.

        Notes
        -----
        This method assumes the image is made by analysis_tools, and
        displays a focal plane plot for LSSTComCam.  Under these
        assumptions, the boxes defined in SVG map map to the detectors
        in the single raft.
        """
        page.extend(['<svg version="1.1" xmlns="http://www.w3.org/2000/svg" '
                     'xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 3600 1800">',
                     f'<image width="3600" height="1800" xlink:href="{relative_file}"></image>',
                     '<a xlink:href="./S00.html"> <rect x="317" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./S01.html"> <rect x="490" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./S02.html"> <rect x="670" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./S10.html"> <rect x="845" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./S11.html"> <rect x="1020" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./S12.html"> <rect x="1197" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./S20.html"> <rect x="1375" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./S21.html"> <rect x="1550" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./S22.html"> <rect x="1550" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '</svg>'])

    def svg_mosaic_latiss(self, relative_file, page):
        """This adds an image to the html, with an SVG overlay for LATISS.

        Parameters
        ----------
        relative_file : `str`
            Relative location of the image to be added.
        page : `list` [`str`]
            The page to append to.

        Notes
        -----
        This method assumes the image is a 8x8 binned focal plane
        mosaic of a LATISS exposure made by one of the VisualizeVisit
        tasks.  Under these assumptions, the boxes defined in SVG map
        to the amplifiers of the single detector.
        """
        page.extend(['<svg version="1.1" xmlns="http://www.w3.org/2000/svg" '
                     'xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 3600 1800">',
                     f'<image width="3600" height="1800" xlink:href="{relative_file}"></image>',
                     '<a xlink:href="./c10.html"> <rect x="317" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c11.html"> <rect x="490" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c12.html"> <rect x="670" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c13.html"> <rect x="845" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c14.html"> <rect x="1020" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c15.html"> <rect x="1197" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c16.html"> <rect x="1375" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c17.html"> <rect x="1550" y="216" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c07.html"> <rect x="1550" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c06.html"> <rect x="1375" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c05.html"> <rect x="1197" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c04.html"> <rect x="1020" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c03.html"> <rect x="845" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c02.html"> <rect x="670" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c01.html"> <rect x="490" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c00.html"> <rect x="317" y="911" fill="#fff" opacity="0" '
                     'width="176" height="692"></rect> </a>',
                     '</svg>'])

    def svg_mosaic_comcam(self, relative_file, page):
        """This adds an image to the html, with an SVG overlay for ComCam.

        Parameters
        ----------
        relative_file : `str`
            Relative location of the image to be added.
        page : `list` [`str`]
            The page to append to.

        Notes
        -----
        This method assumes the image is a 8x8 binned focal plane
        mosaic of a ComCam exposure made by one of the VisualizeVisit
        tasks.  Under these assumptions, the boxes defined in SVG map
        to the detectors in the single raft.
        """
        page.extend(['<svg version="1.1" xmlns="http://www.w3.org/2000/svg" '
                     'xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 640 480">',
                     f'<image width="640" height="480" xlink:href="{relative_file}"></image>',
                     '<a xlink:href="./S00.html"> <rect x="111" y="306" fill="#fff" opacity="0"'
                     'width="128" height="128"></rect> </a>',
                     '<a xlink:href="./S01.html"> <rect x="111" y="175" fill="#fff" opacity="0"'
                     'width="128" height="128"></rect> </a>',
                     '<a xlink:href="./S02.html"> <rect x="111" y="50" fill="#fff" opacity="0" '
                     'width="128" height="128"></rect> </a>',
                     '<a xlink:href="./S10.html"> <rect x="244" y="306" fill="#fff" opacity="0" '
                     'width="128" height="128"></rect> </a>',
                     '<a xlink:href="./S11.html"> <rect x="244" y="175" fill="#fff" opacity="0" '
                     'width="128" height="128"></rect> </a>',
                     '<a xlink:href="./S12.html"> <rect x="244" y="50" fill="#fff" opacity="0" '
                     'width="128" height="128"></rect> </a>',
                     '<a xlink:href="./S20.html"> <rect x="376" y="306" fill="#fff" opacity="0" '
                     'width="128" height="128"></rect> </a>',
                     '<a xlink:href="./S21.html"> <rect x="376" y="175" fill="#fff" opacity="0" '
                     'width="128" height="128"></rect> </a>',
                     '<a xlink:href="./S22.html"> <rect x="376" y="50" fill="#fff" opacity="0" '
                     'width="128" height="128"></rect> </a>',
                     '</svg>'])

    def image_handler(self, image_filename, instrument, page):
        """Handle images.

        This method attempts to find the correct way to add an image
        to a page, with the appropriate SVG subregion information.

        Parameters
        ----------
        image_filename : `str`
            Full path to the image file to be added to a page.
        instrument : `str`
            Instrument name to use to identify special handling.
        page : `list`
            The page to append to.

        """
        relative_file = self.relative_file(image_filename)
        if "Mosaic64" not in relative_file:
            # Skip Mosaic64 products.
            if instrument == 'LATISS':
                if "Mosaic" in relative_file:
                    self.svg_mosaic_latiss(relative_file, page)
                else:
                    self.svg_atools_latiss(relative_file, page)
            elif instrument in ('LSSTComCam', 'LSSTComCamSim'):
                if "Mosaic" in relative_file:
                    # This is a fits_to_png thing
                    self.svg_mosaic_comcam(relative_file, page)
                else:
                    self.svg_atools_comcam(relative_file, page)
            else:
                # We couldn't find a magic handler, so add it
                # manually.
                page.append(f'<img class="center-fit" src="{relative_file}">')

    def block(self, dataset, page, link=None, doc=None):
        """Add a block of information, containing the dataId and collections
        for a particular dataset.

        Parameters
        ----------
        dataset : `dict`
           The dataset containing dataId information.
        page : `list`
           The page to append to.
        link : `str`, optional
           A link to apply to the dataset type name.
        doc : `str`, optional
           A documentation string.  If not supplied, the documentation
           in the configuration is used.
        """
        if doc is None:
            stage = dataset['stage']
            dsType = dataset['type']
            doc = self.dataset_map['stages'][stage][dsType].get("description", "Undocumented.")
        page.extend(["<table width='100%'>",
                     "<tr>",
                     "<th colspan='2' align='center'>DataId:</th>",
                     "<th colspan='2' align='center'>Dataset:</th>", "</tr>"])
        page.extend(["<tr>",
                     "<td align='right'>Instrument</td>", f"<td>{dataset['instrument']}</td>",
                     "<td align='right'>Collection</td>", f"<td>{dataset['collection']}</td>",
                     "</tr>"])
        page.extend(["<tr>",
                     "<td align='right'>Exposure</td>", f"<td>{dataset['exposure']}</td>",
                     "<td align='right'>Stage</td>", f"<td>{dataset['stage']}</td>", "</tr>"])
        page.extend(["<tr>",
                     "<td align='right'>Detector</td>", f"<td>{dataset['detector']}</td>"]),

        page.extend(["<a href='" + f"{link}" + "'>" if link else ""])
        page.extend(["<td align='right'>Dataset type</td>", f"<td>{dataset['type']} ", "</td>"])
        page.extend(["</a>" if link else ""])
        page.extend(["</tr>"])
        page.extend(["<tr>",
                     "<td align='right'>Physical Filter</td>", f"<td>{dataset['physical_filter']}</td>",
                     "<td align='right'>Docstring</td>", f"<td>{doc}</td>",
                     "</tr>"])
        page.extend(["</table>"])

    @staticmethod
    def link(page, text, link):
        text_ref = text.replace(" ", "-")
        page.extend([f"<h2 id='{text_ref}'>", f"<a href='{link}'>",
                     f"{text}" "</a>", "</h2>"])
        return text_ref

    def navigation_init(self):
        page = self.out_files["navigation.html"]
        page.extend(["<table class='ctable' width='100%'>",
                     "<tr>"])

    def navigation_append(self, destination, text, tag, newColumn=False):
        page = self.out_files["navigation.html"]
        if newColumn:
            page.extend(["</tr>", "</tr>"])
        page.extend(["<td>",
                     f"<a href='./{destination}#{tag}'>",
                     f"{text}",
                     "</a>", "</td>"])

    def title(self, page, text):
        text_ref = text.replace(" ", "-")
        page.extend(["<center>",
                     f"<h1 id={text_ref}>", text, "</h1>",
                     "</center>"])
        self.navigation_append("index.html", text, text_ref, newColumn=True)

    def include(self, source, target):
        target_page = self.out_files[target]
        relative_file = self.relative_file(source)
        target_page.extend([f"<iframe width='100%' src='./{relative_file}'>",
                            "</iframe>"])


def main():
    # Could these use the click/daf_butler CLI tools?
    parser = argparse.ArgumentParser(description="Construct a cp_verify metric report.")
    parser.add_argument("-r", "--repository", dest="repository", default="",
                        help="Butler repository to pull results from.")
    parser.add_argument("-O", "--output_path", dest="output_path", default="",
                        help="Output path to write report to.")
    parser.add_argument("-c", "--collections", dest="collections", action="append", default=[],
                        help="Collections to search for results.")
    parser.add_argument("--no_copy", dest="do_copy", action="store_false", default=True,
                        help="Skip copying files for debugging purposes.")
    parser.add_argument("--do_overwrite", dest="do_overwrite", action="store_true", default=False,
                        help="Allow existing files to be overwritten?")
    args = parser.parse_args()

    reporter = CpvReporter(
        repo=args.repository,
        output_path=args.output_path,
        collections=args.collections,
        do_copy=args.do_copy,
        do_overwrite=args.do_overwrite,
    )
    reporter.run()
