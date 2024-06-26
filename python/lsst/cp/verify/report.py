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
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import os
import re

from matplotlib import cm
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from mpl_toolkits.axes_grid1 import make_axes_locatable

from lsst.daf.butler import Butler
from lsst.daf.butler.datastores.file_datastore.retrieve_artifacts import (
    determine_destination_for_retrieved_artifact,
)
from lsst.resources import ResourcePath
from lsst.summit.utils import getQuantiles


class CpvReporter():
    DATASET_MAP = {
        'bias': [
            'bias',
            'biasMosaic8',
            'verifyBiasResults',
            'verifyBiasResidual8',
            'cpBiasCore_metrics',
            'cpBiasCore_biasReadNoisePerAmp_FocalPlaneGeometryPlot',
            'cpBiasCore_biasMeanPerAmp_FocalPlaneGeometryPlot',
            'cpBiasCore_biasNoisePerAmp_FocalPlaneGeometryPlot',
            'cpBiasCore_biasCrNoisePerAmp_FocalPlaneGeometryPlot',
            'cpBiasCore_biasCornerMeanPerAmp_FocalPlaneGeometryPlot',
            'biasIsr_config',
            'verifyBiasApply_config'
        ],
        'dark': [
            'dark',
            'darkMosaic8',
            'verifyDarkResults',
            'verifyDarkResidual8',
            'darkIsr_config',
            'verifyDarkApply_config',
        ],
        'flat': [
            'flat',
            'flatMosaic8',
            'verifyFlatResults',
            'verifyFlatResidual8',
            'flatIsr_config',
            'verifyFlatApply_config',
        ],
        'defects': [
            'defects',
            'verifyDefectsResults',
        ],
        'ptc': [
            'ptc',
            'verifyPtcResults',
        ],
        'linearity': [],
        'crosstalk': [],
        'bfk': [],
        'cti': [],
    }

    DOC_MAP = {
        'bias': "The combined bias calibration.",
        'biasMosaic8': "Mosaic of combined bias calibration.",
        'verifyBiasResults': "Catalog of combined bias verification results.",
        'verifyBiasResidual8': "Mosaic of bias residuals (bias exposure corrected up through bias application).",
        'cpBiasCore_metrics': "Core metric bundle for bias calibration.",
        'cpBiasCore_biasCornerMeanPerAmp_FocalPlaneGeometryPlot': "The bias mean in 200x200 box at amp readout corner.",
        'cpBiasCore_biasCrNoisePerAmp_FocalPlaneGeometryPlot': "The image noise measured after cosmic ray rejection.",
        'cpBiasCore_biasMeanPerAmp_FocalPlaneGeometryPlot': "The image mean.",
        'cpBiasCore_biasNoisePerAmp_FocalPlaneGeometryPlot': "The image noise.",
        'cpBiasCore_biasReadNoisePerAmp_FocalPlaneGeometryPlot': "The noise in the serial overscan after overscan correction.",
        'dark': "The combined dark calibration.",
        'darkMosaic8': "Mosaic of combined dark calibration.",
        'verifyDarkResults': "Catalog of combined dark verification results.",
        'verifyDarkResidual8': "Mosaic of dark residuals (dark exposure corrected up through dark application).",
        'flat': "The combined flat calibration.",
        'flatMosaic8': "Mosaic of combined flat calibration.",
        'verifyFlatResults': "Catalog of combined flat verification results.",
        'verifyFlatResidual8': "Mosaic of flat residuals (flat exposure corrected up through flat application).",
        'defects': "The combined defects calibration.",
        'verifyDefectsResults': "Catalog of combined defect verification results.",
        'ptc': "The photon transfer curve calibration.",
        'verifyPtcResults': "Catalog of PTC verification results."
    }

    if False:
        DESTINATION = "/sdf/home/c/czw/public_html/cpv_reports/PREOPS-5181/"
        COLLECTIONS = [
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyBias.20240529a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyDark.20240529a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyFlat-r.20240529a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyFlat-g.20240529a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyFlat-z.20240529a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyFlat-y.20240529a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyFlat-white.20240529a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyBias.20240528a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyDark.20240528a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyFlat-r.20240528a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyFlat-g.20240528a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyFlat-z.20240528a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyFlat-y.20240528a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyFlat-empty.20240528a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/ptcGen.20240528a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyDark.20240530a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyFlat-r.20240530a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyFlat-g.20240530a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyFlat-z.20240530a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyFlat-y.20240530a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/verifyFlat-white.20240530a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/ptcGen.20240530a",
            "u/huanlin/PREOPS-5181/CR1_LATISS/ptcGen.20240530a.ABC",
        ]
    else:
        DESTINATION = "/sdf/home/c/czw/public_html/cpv_reports/PREOPS-5182/"
        COLLECTIONS = [
            "u/huanlin/PREOPS-5182/CR1_LSSTComCamSim/verifyBias.20240528a",
            "u/huanlin/PREOPS-5182/CR1_LSSTComCamSim/verifyDark.20240528a",
            "u/huanlin/PREOPS-5182/CR1_LSSTComCamSim/verifyFlat-i06.20240528a",
            "u/huanlin/PREOPS-5182/CR1_LSSTComCamSim/verifyBias.20240530a",
            "u/huanlin/PREOPS-5182/CR1_LSSTComCamSim/verifyDark.20240530a",
            "u/huanlin/PREOPS-5182/CR1_LSSTComCamSim/verifyFlat-i06.20240530a",
        ]
    DO_COPY = True
    DO_OVERWRITE = False
    def __init__(self, repo="/repo/embargo", **kwargs):
        super().__init__()
        self.datasets = []
        self.out_files = {
            'index.html': self._init_page(),
        }
        self.repo = repo
        self.butler = Butler(repo)

    def generate_report(self):
        src_dir = os.path.join(self.DESTINATION, "src")
        os.makedirs(src_dir, exist_ok=True)

        self.copy_datasets(doCopy=self.DO_COPY)
        self.parse_datasets()
        self.write_pages()

    def copy_datasets(self, doCopy=None, doOverwrite=None,
                      collections=None,
                      ):
        doCopy = self.DO_COPY
        doOverwrite = self.DO_OVERWRITE
        collections = self.COLLECTIONS
        # Based on https://github.com/lsst/daf_butler/blob/main/python/lsst/daf/butler/datastores/fileDatastore.py#L1978 # noqa W505
        for stage, dataset_types in self.DATASET_MAP.items():
            refs = self.butler.registry.queryDatasets(datasetType=dataset_types,
                                                      collections=collections,
                                                      where=None, findFirst=True)
            for ref in refs:
                locations = self.butler._datastore._get_dataset_locations_info(ref)
                print(f"Found {stage} {ref}")
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
                for location, _ in locations:
                    source_uri = location.uri
                    target_uri = determine_destination_for_retrieved_artifact(
                        ResourcePath(os.path.join(self.DESTINATION, "src")),
                        location.pathInStore,
                        False
                    )
                    dataset['uri'] = target_uri
                    willCopy = doCopy

                    if target_uri.getExtension().lower() == ".fits" and \
                       dataset['storage_class'] in ('ImageF', 'ExposureF', 'ImageD', 'ExposureD'):
                        # We need a PNG at that location.
                        png_uri = target_uri.updatedExtension("png")
                        dataset['uri'] = png_uri
                        if os.path.exists(png_uri.path):
                            if willCopy and not doOverwrite:
                                willCopy = False

                        if willCopy:
                            data = self.butler.get(ref)
                            self.fits_to_png(data, png_uri.path, f"{ref.datasetType.name} {ref.dataId}")
                    elif target_uri.getExtension().lower() == ".parq" and \
                         dataset['storage_class'] in ('ArrowAstropy', ):
                        # Let's convert the table for now
                        html_uri = target_uri.updatedExtension("html")
                        dataset['uri'] = html_uri
                        if os.path.exists(html_uri.path):
                            if willCopy and not doOverwrite:
                                willCopy = False
                        willCopy=True
                        if willCopy:
                            data = self.butler.get(ref)
                            files_created = self.parq_to_html(data, html_uri)
                    elif target_uri.getExtension().lower() == ".py":
                        # These should be config files
                        html_uri = target_uri.updatedExtension("html")
                        dataset['uri'] = html_uri
                        if willCopy:
                            data = self.butler.get(ref)
                            self.py_to_html(data, html_uri)
                    elif os.path.exists(target_uri.path):
                        if willCopy and doOverwrite:
                            # Otherwise just copy directly.
                            target_uri.transfer_from(source_uri, transfer="copy", overwrite=True)
                    else:
                        if willCopy:
                            # Otherwise just copy directly.
                            target_uri.transfer_from(source_uri, transfer="copy", overwrite=True)

                self.datasets.append(dataset)

    def parse_datasets(self):
        populated_stages_list = []

        # Make navigation page.
        self.out_files["navigation.html"] = self._init_page("Navigation")
        self.navigation_init()
        self.include("navigation.html", "index.html")
        self.out_files["manifest.html"] = self._init_page("Manifest")

        tag = self.link(self.out_files["index.html"], "Parent", "./..")
        self.navigation_append("index.html", "Parent", tag, newColumn=True)
        self.navigation_append("manifest.html", "Manifest", "", newColumn=False)

        for stage in ['bias', 'dark', 'flat', 'ptc', 'defects', 'crosstalk', 'cti']:
            dss = [x for x in self.datasets if x['stage'] == stage]
            if len(dss) > 0:
                populated_stages_list.append(stage)
            else:
                continue
            self.title(self.out_files["index.html"], stage)

            self.out_files[f"{stage}_exp.html"] = self._init_page()
            self.include("navigation.html", f"{stage}_exp.html")
            self.out_files[f"{stage}_det.html"] = self._init_page()
            self.include("navigation.html", f"{stage}_det.html")

            exps = [x for x in dss if x['exposure'] is not None]
            dets = [x for x in dss if x['detector'] is not None]

            for ds in sorted(exps, key=lambda x: x['exposure']):
                self.block(ds, self.out_files[f"{stage}_exp.html"])
                if ('.png' in ds['uri'].path):
                    self.image_handler(ds['uri'].path, self.out_files[f"{stage}_exp.html"])
            tag = self.link(self.out_files["index.html"], f"{stage} exposures", f"./{stage}_exp.html")
            self.navigation_append("index.html", f"{stage} exposures", tag)

            for ds in sorted(dets, key=lambda x: x['detector']):
                self.block(ds, self.out_files[f"{stage}_det.html"])
                if ('.png' in ds['uri'].path):
                    self.image_handler(ds['uri'].path, self.out_files[f"{stage}_det.html"])
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
                        self.image_handler(ds['uri'].path, self.out_files["index.html"])
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


    def write_pages(self):
        for ff, contents in self.out_files.items():
            reportOut = os.path.join(self.DESTINATION, ff)
            contents.extend(self._close_page())
            with open(reportOut, 'w') as ff:
                for line in contents:
                    print(line, file=ff)

    @staticmethod
    def _init_page(title=""):
        return [
            "<html>", "<head>", "<style> ",
            " * { margin: 0; padding: 0;}",
            " .imgbox { display: grid; height: 100%; }",
            " .center-fit { max-width: 100%; max-height: 100vh; margin:auto; }",
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
            '<map name="atools_comcam">',
            '<area shape="rect" coords="325,1600,780,1155" alt="S00" href="S00.html">',
            '<area shape="rect" coords="325,1130,780,685" alt="S01" href="S01.html">',
            '<area shape="rect" coords="325,660,780,215" alt="S02" href="S02.html">',
            '<area shape="rect" coords="795,1600,1250,1155" alt="S00" href="S10.html">',
            '<area shape="rect" coords="795,1130,1250,685" alt="S01" href="S11.html">',
            '<area shape="rect" coords="795,660,1250,215" alt="S02" href="S12.html">',
            '<area shape="rect" coords="1265,1600,1720,1155" alt="S00" href="S20.html">',
            '<area shape="rect" coords="1265,1130,1720,685" alt="S01" href="S21.html">',
            '<area shape="rect" coords="1265,660,1720,215" alt="S02' href="S22.html'>',
            '</map>',
            '<map name="atools_latiss">',
            '<area shape="rect" coords="
 <!-- (317,1600) @ 490, 670, 845, 1020, 1197, 1375, 1550, 1725 -->
 <!--            @ 909, 1380 -->
</map>
<map name="mosaic_latiss">
 <!-- (115, 431) - (500, 50) // 8, 2 -->
</map>
<map name="mosaic_comcam">
 <!-- (111, 437) - (240, 313) , (244, 437) - (372, 314), (376, 438) - (503, 313) -->
 <!--                                                          306) - (504, 180) -->
 <!-- 							       175) - (504, 47)  -->
</map>
        ]

    @staticmethod
    def _close_page():
        return ["</body>",
                "</html>"]

    @staticmethod
    def relative_file(filename):
        return re.sub(r"^.*/src/", "./src/", filename)

    def image_handler(self, image_filename, page):
        """Handle images."""
        relative_file = self.relative_file(image_filename)
        if "Mosaic64" not in relative_file:
            page.append(f'<a href="{relative_file}">')
            page.append(f'<img class="center-fit" src="{relative_file}">')
            page.append('</a>')


    def block(self, dataset, page, link=None, doc=None):
        if doc is None:
            doc = self.DOC_MAP.get(dataset['type'], "Undocumented.")
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
        if link:
            page.extend(["<a href='" + f"{link}" + "'>",])
        page.extend(["<td align='right'>Dataset type</td>", f"<td>{dataset['type']} ", "</td>"])
        if link:
            page.extend(["</a>"])
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

    @staticmethod
    def fits_to_png(data, target_uri, title):
        # Taken from
        # https://github.com/lsst-sitcom/rubintv_production/blob/main/python/lsst/rubintv/production/slac/mosaicing.py  # noqa W505
        # _plotFpMosaic()
        try:
            array = data.image.array
        except AttributeError:
            array = data.array

        # This should be available in lsst.utils.plotting, as of DM-44725.
        fig =  Figure()
        canvas = FigureCanvasAgg(fig)
        ax = fig.gca()
        ax.clear()
        cmap = cm.gray
        quantiles = getQuantiles(array, cmap.N)
        norm = colors.BoundaryNorm(quantiles, cmap.N)

        im = ax.imshow(array, norm=norm, interpolation='None', cmap=cmap, origin='lower')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.suptitle(title)
        fig.colorbar(im, cax=cax)
        fig.tight_layout()
        fig.savefig(target_uri.path)
        plt.close()

    @staticmethod
    def parq_to_html(data, target_uri):
        format_dict = {}

        files_created = {}
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
        # files_created[f"{filename_base}.html"] = {}

        if False:
            # Write subset tables:
            if 'exposure' in data.columns:
                exposures = set(data['exposure'])
                for exp in exposures:
                    mask = data['exposure'] == exp
                    subset = data[mask]
                    data.write(f"{filename_base}_exp{exp}.html",
                               format='ascii.html',
                               overwrite=True,
                    formats=format_dict)
                    files_created[f"{filename_base}_exp{exp}.html"] = {'exposure': exp}

                    if 'detector' in data.columns:
                        detectors = set(data['exposure'])
                        for det in detectors:
                            mask = data['detector'] == det
                            subset = data[mask]
                            data.write(f"{filename_base}_det{det}.html",
                                       format='ascii.html',
                                       overwrite=True,
                                       formats=format_dict)
                            files_created[f"{filename_base}_det{det}.html"] = {'detector': det}

    def py_to_html(self, data, target_uri):
        # For configs and such.
        with open(target_uri.path, 'w') as ff:
            for line in self._init_page():
                print(line, file=ff)
                # import pdb; pdb.set_trace()
            for line in data.saveToString().split("\n"):
                if len(line) > 0 and  line[0] == '#':
                    print("<pre class='pre-comment'>", line, "</pre>", file=ff)
                else:
                    print("<pre>", line, "</pre>", file=ff)
            for line in self._close_page():
                print(line, file=ff)


def main():
    reporter = CpvReporter()
    reporter.generate_report()
