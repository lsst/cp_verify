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
from mpl_toolkits.axes_grid1 import make_axes_locatable

# import lsst.afw.image as afwImage

from lsst.daf.butler import Butler
from lsst.daf.butler.datastores.file_datastore.retrieve_artifacts import (
    determine_destination_for_retrieved_artifact,
)
from lsst.resources import ResourcePath
from lsst.summit.utils import getQuantiles


class CpvReporter():
    DESTINATION = "/sdf/home/c/czw/public_html/cpv_reports/PREOPS-5182/"

    DATASET_MAP = {
        'bias': [
            'bias',
            'biasMosaic8',
            'biasMosaic64',
            'verifyBiasResults',
            'verifyBiasResidual8',
            'verifyBiasResidual64',
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
            'darkMosaic64',
            'verifyDarkResults',
            'verifyDarkResidual8',
            'verifyDarkResidual64',
        ],
        'flat': [
            'flat',
            'flatMosaic8',
            'flatMosaic64',
            'verifyFlatResults',
            'verifyFlatResidual8',
            'verifyFlatResidual64',
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

    def __init__(self, repo="/repo/embargo", **kwargs):
        super().__init__()
        self.datasets = []
        self.out_files = {
            'index.html': self._init_page(),
        }
        self.repo = repo
        self.butler = Butler(repo)

    def copy_datasets(self, doCopy=True,
                      collections=[
                          "u/huanlin/PREOPS-5182/CR1_LSSTComCamSim/verifyBias.20240528a",
                          "u/huanlin/PREOPS-5182/CR1_LSSTComCamSim/verifyDark.20240528a",
                          "u/huanlin/PREOPS-5182/CR1_LSSTComCamSim/verifyFlat-i06.20240528a",
                          #                          "u/huanlin/PREOPS-5181/CR1_LATISS/verifyBias.20240529a",

                      ]):
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

                    if target_uri.getExtension().lower() == ".fits":
                        # We need a PNG at that location.
                        png_uri = target_uri.updatedExtension("png")
                        dataset['uri'] = png_uri

                        if doCopy:
                            data = self.butler.get(ref)
                            self.fits_to_png(data, png_uri.path, f"{ref.datasetType.name} {ref.dataId}")
                    elif target_uri.getExtension().lower() == ".parq":
                        # Let's convert the table for now
                        html_uri = target_uri.updatedExtension("html")
                        dataset['uri'] = html_uri
                        if doCopy:
                            data = self.butler.get(ref)
                            data.write(html_uri.path, format='ascii.html')
                    elif doCopy:
                        # Otherwise just copy directly.
                        target_uri.transfer_from(source_uri, transfer="copy", overwrite=True)

                self.datasets.append(dataset)

    def parse_datasets(self):
        populated_stages_list = []

        for stage in ['bias', 'dark', 'flat', 'ptc', ]:
            dss = [x for x in self.datasets if x['stage'] == stage]
            if len(dss) > 0:
                populated_stages_list.append(stage)

            self.out_files[f"{stage}.html"] = self._init_page()
            self.out_files[f"{stage}_exp.html"] = self._init_page()
            self.out_files[f"{stage}_det.html"] = self._init_page()

            exps = [x for x in dss if x['exposure'] is not None]
            dets = [x for x in dss if x['detector'] is not None]

            for ds in sorted(exps, key=lambda x: x['exposure']):
                if ('.png' in ds['uri'].path):
                    self.block(ds, self.out_files[f"{stage}_exp.html"])
                    self.image_handler(ds['uri'].path, self.out_files[f"{stage}_exp.html"])
            self.link(self.out_files["index.html"], f"{stage} exposures", f"./{stage}_exp.html")

            for ds in sorted(dets, key=lambda x: x['detector']):
                if ('.png' in ds['uri'].path):
                    self.block(ds, self.out_files[f"{stage}_det.html"])
                    self.image_handler(ds['uri'].path, self.out_files[f"{stage}_det.html"])
            self.link(self.out_files["index.html"], f"{stage} detectors", f"./{stage}_det.html")

            for ds in dss:
                if ds['storageClass'] == 'ArrowAstropy':
                    # Is a results catalog
                    self.block(ds, self.out_files["index.html"], ds['uri'])
                elif ds['exposure'] is None and ds['detector'] is None:
                    self.block(ds, self.out_files["index.html"])
                    self.image_handler(ds['uri'].path, self.out_files["index.html"])
            self.link(self.out_files["index.html"], f"{stage}", f"./{stage}.html")

    def generate_report(self):
        src_dir = os.path.join(self.DESTINATION, "src")
        os.makedirs(src_dir, exist_ok=True)

        self.copy_datasets(doCopy=True)
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
            "</style>",
            "</title>", title, "</title>",
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
        page.append(f'<a href="{relative_file}">')
        page.append(f'<img class="center-fit" src="{relative_file}">')
        page.append('</a>')

    @staticmethod
    def block(dataset, page, link=None):
        page.extend(["<table width='100%'>", "<tr>", "<th>DataId</th>", "<th>T1</th>", "</tr>"])
        page.extend(["<tr>",
                     f"<td>{dataset['instrument']}</td>",
                     f"<td>{dataset['collection']}</td>",
                     "</tr>"])
        page.extend(["<tr>", f"<td>{dataset['exposure']}</td>", f"<td>{dataset['stage']}</td>", "</tr>"])
        if link:
            page.extend(["<tr>", f"<td>{dataset['detector']}</td>"]),
            page.extend(["<a href='" + f"{dataset['uri']}" + "'>",
                         "<td>{dataset['type']}</td>", "</a>", "</tr>"])
        else:
            page.extend(["<tr>", f"<td>{dataset['detector']}</td>", f"<td>{dataset['type']}</td>", "</tr>"])
        page.extend(["<tr>", f"<td>{dataset['physical_filter']}</td>", "<td></td>", "</tr>"])
        page.extend(["</table>"])

    @staticmethod
    def link(page, text, link):
        page.extend(["<table width='100%'>", "<tr>", "<th>T0</th>", "<th>T1</th>", "</tr>"])
        page.extend(["<tr>", "<td></td>"]),
        page.extend(["<td>", f"<a href='{link}'>",
                     f"{text}" "</a>", "</td>", "</tr>"])
        page.extend(["</table>"])

    @staticmethod
    def fits_to_png(data, out_filename, title):
        # Taken from
        # https://github.com/lsst-sitcom/rubintv_production/blob/main/python/lsst/rubintv/production/slac/mosaicing.py  # noqa W505
        # _plotFpMosaic()
        try:
            array = data.image.array
        except AttributeError:
            array = data.array

        fig = plt.figure()
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
