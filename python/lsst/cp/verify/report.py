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

# import glob
# import json
# import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import os
import re

from matplotlib import cm
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from mpl_toolkits.axes_grid1 import make_axes_locatable

# import lsst.afw.image as afwImage

from lsst.daf.butler import Butler
from lsst.daf.butler.datastores.file_datastore.retrieve_artifacts import (
    determine_destination_for_retrieved_artifact,
)
from lsst.resources import ResourcePath
from lsst.summit.utils import getQuantiles
# from fgcm.fgcmUtilities import makeFigure  # This doesn't find it??

class CpvReporter():
    DATASET_MAP = {
        'bias': [
            'bias',
            'biasMosaic8',
            'verifyBiasResults',
            'verifyBiasResidual8',
            'cpBiasCore_metrics',
            'cpBiasCore_biasCornerMeanPerAmp_FocalPlaneGeometryPlot',
            'cpBiasCore_biasCrNoisePerAmp_FocalPlaneGeometryPlot',
            'cpBiasCore_biasMeanPerAmp_FocalPlaneGeometryPlot',
            'cpBiasCore_biasNoisePerAmp_FocalPlaneGeometryPlot',
            'cpBiasCore_biasReadNoisePerAmp_FocalPlaneGeometryPlot',
        ],
        'dark': [
            'dark',
            'darkMosaic8',
            #            'darkMosaic64',
            'verifyDarkResults',
            'verifyDarkResidual8',
            # 'verifyDarkResidual64',
        ],
        'flat': [
            'flat',
            'flatMosaic8',
            # 'flatMosaic64',
            'verifyFlatResults',
            'verifyFlatResidual8',
            # 'verifyFlatResidual64',
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


    if True:
        DESTINATION = "/sdf/home/c/czw/public_html/cpv_reports/PREOPS-5181/"
        DESTINATION = "/sdf/home/c/czw/public_html/cpv_reports/PREOPS-5181a/"
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
        DESTINATION = "/sdf/home/c/czw/public_html/cpv_reports/2024-06-06/"
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
                           'storageClass': ref.datasetType.storageClass.name,
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

                    if target_uri.getExtension().lower() == ".fits" and dataset['storageClass'] in ('ImageF', 'ExposureF', 'ImageD', 'ExposureD'):
                        # We need a PNG at that location.
                        png_uri = target_uri.updatedExtension("png")
                        dataset['uri'] = png_uri
                        if os.path.exists(png_uri.path):
                            if willCopy and not doOverwrite:
                                willCopy = False

                        if willCopy:
                            data = self.butler.get(ref)
                            self.fits_to_png(data, png_uri.path, f"{ref.datasetType.name} {ref.dataId}")
                    elif target_uri.getExtension().lower() == ".parq":
                        # Let's convert the table for now
                        html_uri = target_uri.updatedExtension("html")
                        dataset['uri'] = html_uri
                        if os.path.exists(html_uri.path):
                            if willCopy and not doOverwrite:
                                willCopy = False
                        willCopy=True
                        if willCopy:
                            data = self.butler.get(ref)
                            data.write(html_uri.path, format='ascii.html', overwrite=True,
                                       exclude_names=['BIAS_SERIAL_PROF', 'BIAS_PARALLEL_PROF'] )

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
        self.include("navigation.html", "index.html")
        self.out_files["manifest.html"] = self._init_page("Manifest")

        tag = self.link(self.out_files["index.html"], "Parent", "./..")
        self.navigation_append("index.html", "Parent", tag)

        for stage in ['bias', 'dark', 'flat', 'ptc', ]:
            dss = [x for x in self.datasets if x['stage'] == stage]
            if len(dss) > 0:
                populated_stages_list.append(stage)
            self.title(self.out_files["index.html"], stage)

            # self.out_files[f"{stage}.html"] = self._init_page()
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
                if ds['storageClass'] == 'ArrowAstropy':
                    # Is a results catalog
                    self.block(ds, self.out_files["index.html"], ds['uri'])
                elif ds['exposure'] is None and ds['detector'] is None:
                    self.block(ds, self.out_files["index.html"])
                    if ".png" in ds['uri'].path:
                        self.image_handler(ds['uri'].path, self.out_files["index.html"])

            #self.link(self.out_files["index.html"], f"{stage}", f"./{stage}.html")

        # Record all datasets here, and link to them directly:
        page = self.out_files["manifest.html"]
        tag = self.link(self.out_files["index.html"], "Manifest", f"./manifest.html")
        self.navigation_append("index.html", "Manifest", tag)

        page.extend(["<table width='100%'>",
                     "<tr>",
                     "<th>stage</th>", "<th>ref</th>", "<th>type</th>", "<th>storageClass</th>",
                     "<th>dataId</th>", "<th>collection</th>", 
                     "</tr>"])
        for ds in self.datasets:
            relative_file = re.sub(r"^.*/src/", "./src/", ds['uri'].path)
            page.extend(["<tr>",
                         f"<td>{ds['stage']}</td>", f"<td>{ds['ref']}</td>",
                         "<td>", f'<a href="{relative_file}">', f"{ds['type']}", "</a>",
                         f"<td>{ds['storageClass']}</td>",
                         f"<td>{ds['dataId']}</td>", f"<td>{ds['collection']}</td>",
                         "</tr>"])


    def generate_report(self):
        src_dir = os.path.join(self.DESTINATION, "src")
        os.makedirs(src_dir, exist_ok=True)

        self.copy_datasets(doCopy=self.DO_COPY)
        self.parse_datasets()
        self.write_pages()

    def write_pages(self):
        for ff, contents in self.out_files.items():
            reportOut = os.path.join(self.DESTINATION, ff)
            contents.extend(self._close_page())
            with open(reportOut, 'w') as ff:
                for line in contents:
                    print(line, file=ff)

    @staticmethod
    def _init_page(title=None):
        return [
            "<html>", "<head>", "<style> ",
            " * { margin: 0; padding: 0;}",
            " .imgbox { display: grid; height: 100%; }",
            " .center-fit { max-width: 100%; max-height: 100vh; margin:auto; }",
            " .tooltip { position: relative; display: inline-block; border-bottom: 1px dotted black; }",
            " .tooltip .tooltiptext { visibility: hidden; color: white; background-color: black; text-align: center; top: -5px; left: 105%; }",
            " .tooltip:hover .tooltiptext { visibility: visible; }",
            "</style>",
            "</title>", title, "</title>",
            '<base target="_parent">',
            "</head>",
            "<body>",
        ]

    @staticmethod
    def _close_page():
        return ["</body>",
                "</html>"]

    @staticmethod
    def image_handler(image_filename, page):
        """Handle images."""
        relative_file = re.sub(r"^.*/src/", "./src/", image_filename)
        if "Mosaic64" not in relative_file:
            page.append(f'<a href="{relative_file}">')
            page.append(f'<img class="center-fit" src="{relative_file}">')
            page.append('</a>')


    def block(self, dataset, page, link=None, doc=None):
        if doc is None:
            doc = self.DOC_MAP.get(dataset['type'], "Undocumented.")
        page.extend(["<table width='100%'>",
                     "<tr>",
                     "<th></th>", "<th align='left'>DataId:</th>",
                     "<th></th>", "<th align='left'>Dataset:</th>", "</tr>"])
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

    def navigation_append(self, destination, text, tag, level=3):
        page = self.out_files["navigation.html"]
        page.extend([f"<a href='./{destination}#{tag}'>",
                     f"<h{level}>{text}</h{level}>",
                     "</a>"])

    def title(self, page, text):
        text_ref = text.replace(" ", "-")
        page.extend(["<center>",
                     f"<h1 id={text_ref}>", text, "</h1>",
                     "</center>"])
        self.navigation_append("index.html", text, text_ref)

    def include(self, source, target):
        source_page = self.out_files[source]
        target_page = self.out_files[target]

        target_page.extend([f"<iframe width='100%' src='./{source}'>",
                            "</iframe>"])

    @staticmethod
    def fits_to_png(data, out_filename, title):
        # Taken from
        # https://github.com/lsst-sitcom/rubintv_production/blob/main/python/lsst/rubintv/production/slac/mosaicing.py  # noqa W505
        # _plotFpMosaic()
        try:
            array = data.image.array
        except AttributeError:
            array = data.array

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
        fig.savefig(out_filename)
        plt.close()

    @staticmethod
    def table_handler(data_filename, page):
        """Handle data tables."""
        pass


def main():
    reporter = CpvReporter()
    reporter.generate_report()
