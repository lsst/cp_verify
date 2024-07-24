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
            'cpPtcCore_metrics',
            'cpPtcCore_ptcGainPerAmp_FocalPlaneGeometryPlot',
            'cpPtcCore_ptcNoisePerAmp_FocalPlaneGeometryPlot',
            'cpPtcCore_ptcA00PerAmp_FocalPlaneGeometryPlot',
            'cpPtcCore_ptcTurnoffPerAmp_FocalPlaneGeometryPlot',
            'cpPtcCore_ptcTurnoffPerAmp_FocalPlaneGeometryPlot',
            'cpPtcCore_ptcRowMeanVarianceSlopePerAmp_FocalPlaneGeometryPlot',
            'cpPtcDetCore_ptcPlot_GridPlot',

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
        'verifyPtcResults': "Catalog of PTC verification results.",
        'cpPtcCore_metrics': "Core metric bundle for ptc calibration.",
        'cpPtcCore_ptcGainPerAmp_FocalPlaneGeometryPlot': "TBA",
        'cpPtcCore_ptcNoisePerAmp_FocalPlaneGeometryPlot': "TBA",
        'cpPtcCore_ptcA00PerAmp_FocalPlaneGeometryPlot': "TBA",
        'cpPtcCore_ptcTurnoffPerAmp_FocalPlaneGeometryPlot': "TBA",
        'cpPtcCore_ptcTurnoffPerAmp_FocalPlaneGeometryPlot': "TBA",
        'cpPtcCore_ptcRowMeanVarianceSlopePerAmp_FocalPlaneGeometryPlot': "TBA",
        'cpPtcDetCore_ptcPlot': "TBA",
        'cpPtcDetCore_ptcPlot_GridPlot': "TBA",
    }

    REPO = "/repo/embargo"
    OPRE = 4
    INST = 'LSSTComCamSim'
    # INST = 'LATISS'
    if (OPRE == 4) and (INST == 'LATISS'):
        DESTINATION = "/sdf/home/c/czw/public_html/cpv_reports/PREOPS-5261"
        COLLECTIONS = [
            "u/huanlin/PREOPS-5261/OR4_LATISS/verifyBias.20240625a",
            "u/huanlin/PREOPS-5261/OR4_LATISS/verifyDark.20240625a",
            "u/huanlin/PREOPS-5261/OR4_LATISS/verifyFlat-r.20240625a",
            "u/huanlin/PREOPS-5261/OR4_LATISS/verifyFlat-g.20240625a",
            "u/huanlin/PREOPS-5261/OR4_LATISS/verifyFlat-z.20240625a",
            "u/huanlin/PREOPS-5261/OR4_LATISS/verifyFlat-y.20240625a",
            "u/huanlin/PREOPS-5261/OR4_LATISS/verifyFlat-empty.20240625a",
            "u/huanlin/PREOPS-5261/OR4_LATISS/ptcGen.20240625a.A",
            "u/huanlin/PREOPS-5261/OR4_LATISS/verifyBias.20240626a",
            "u/huanlin/PREOPS-5261/OR4_LATISS/verifyDark.20240626a",
            "u/huanlin/PREOPS-5261/OR4_LATISS/verifyFlat-r.20240626a",
            "u/huanlin/PREOPS-5261/OR4_LATISS/verifyFlat-g.20240626a",
            "u/huanlin/PREOPS-5261/OR4_LATISS/verifyFlat-z.20240626a",
            "u/huanlin/PREOPS-5261/OR4_LATISS/verifyFlat-y.20240626a",
            "u/huanlin/PREOPS-5261/OR4_LATISS/verifyFlat-white.20240626a",
            "u/huanlin/PREOPS-5261/OR4_LATISS/ptcGen.20240626a.B",
        ]
    elif (OPRE == 4) and (INST == 'LSSTComCamSim'):
        REPO = 'embargo_or4'
        DESTINATION = "/sdf/home/c/czw/public_html/cpv_reports/PREOPS-5262"
        COLLECTIONS = [
            "u/huanlin/PREOPS-5262/OR4_LSSTComCamSim/verifyBias.20240625a",
            "u/huanlin/PREOPS-5262/OR4_LSSTComCamSim/verifyDark.20240625a",
            "u/huanlin/PREOPS-5262/OR4_LSSTComCamSim/verifyFlat-i06.20240625a",
        ]
    elif (OPRE == 'cal') and (INST == 'LATISS'):
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
    elif (OPRE == 'cal') and (INST == 'LSSTComCamSim'):
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
    
    def __init__(self, repo=REPO, **kwargs):
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
                            self.fits_to_png(data, png_uri, f"{ref.datasetType.name} {ref.dataId}")
                    elif target_uri.getExtension().lower() == ".parq" and dataset['storage_class'] in ('ArrowAstropy', ):
                        # Let's convert the table for now
                        html_uri = target_uri.updatedExtension("html")
                        dataset['uri'] = html_uri
                        if os.path.exists(html_uri.path):
                            if willCopy and not doOverwrite:
                                willCopy = False
                        willCopy = True
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
            # " .center-fit { max-width: 100%; max-height: 100vh; margin:auto; }",
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
            '<area shape="rect" coords="1265,660,1720,215" alt="S02" href="S22.html">',
            '</map>',

            '<map name="atools_latiss">',
            '<area shape="rect" coords="317,1400,490,1380" alt="C10" href="C10.html">',
            '<area shape="rect" coords="490,1600,670,1380" alt="C11" href="C11.html">',
            '<area shape="rect" coords="670,1600,845,1380" alt="C12" href="C12.html">',
            '<area shape="rect" coords="845,1600,1020,1380" alt="C13" href="C13.html">',
            '<area shape="rect" coords="1020,1600,1197,1380" alt="C14" href="C14.html">',
            '<area shape="rect" coords="1197,1600,1376,1380" alt="C15" href="C15.html">',
            '<area shape="rect" coords="1375,1600,1550,1380" alt="C16" href="C16.html">',
            '<area shape="rect" coords="1550,1600,1725,1380" alt="C17" href="C17.html">',
            '<area shape="rect" coords="317,1380,490,909" alt="C00" href="C00.html">',
            '<area shape="rect" coords="490,1380,670,909" alt="C01" href="C01.html">',
            '<area shape="rect" coords="670,1380,845,909" alt="C02" href="C02.html">',
            '<area shape="rect" coords="845,1380,1020,909" alt="C03" href="C03.html">',
            '<area shape="rect" coords="1020,1380,1197,909" alt="C04" href="C04.html">',
            '<area shape="rect" coords="1197,1380,1376,909" alt="C05" href="C05.html">',
            '<area shape="rect" coords="1375,1380,1550,909" alt="C06" href="C06.html">',
            '<area shape="rect" coords="1550,1380,1725,909" alt="C07" href="C07.html">',
            '</map>'

            '<map name="mosaic_latiss">',
            '<area shape="rect" coords="115,431,163,241" alt="C10" href="C10.html">',
            '<area shape="rect" coords="163,431,211,241" alt="C11" href="C11.html">',
            '<area shape="rect" coords="211,431,260,241" alt="C12" href="C12.html">',
            '<area shape="rect" coords="260,431,307,241" alt="C13" href="C13.html">',
            '<area shape="rect" coords="307,431,356,241" alt="C14" href="C14.html">',
            '<area shape="rect" coords="356,431,404,241" alt="C15" href="C15.html">',
            '<area shape="rect" coords="404,431,451,241" alt="C16" href="C16.html">',
            '<area shape="rect" coords="451,431,500,241" alt="C17" href="C17.html">',
            '<area shape="rect" coords="115,241,163,50" alt="C00" href="C00.html">',
            '<area shape="rect" coords="163,241,211,50" alt="C01" href="C01.html">',
            '<area shape="rect" coords="211,241,260,50" alt="C02" href="C02.html">',
            '<area shape="rect" coords="260,241,307,50" alt="C03" href="C03.html">',
            '<area shape="rect" coords="307,241,356,50" alt="C04" href="C04.html">',
            '<area shape="rect" coords="356,241,404,50" alt="C05" href="C05.html">',
            '<area shape="rect" coords="404,241,451,50" alt="C06" href="C06.html">',
            '<area shape="rect" coords="451,241,500,50" alt="C07" href="C07.html">',
            '</map>'

            '<map name="mosaic_comcam">',
            '<area shape="rect" coords="111,175,240,47" alt="S00" href="S00.html">',
            '<area shape="rect" coords="111,175,240,47" alt="S01" href="S01.html">',
            '<area shape="rect" coords="111,175,240,47" alt="S02" href="S02.html">',
            '<area shape="rect" coords="244,306,372,180" alt="S00" href="S10.html">',
            '<area shape="rect" coords="244,306,372,180" alt="S01" href="S11.html">',
            '<area shape="rect" coords="244,306,372,180" alt="S02" href="S12.html">',
            '<area shape="rect" coords="376,438,503,313" alt="S00" href="S20.html">',
            '<area shape="rect" coords="376,438,503,313" alt="S01" href="S21.html">',
            '<area shape="rect" coords="376,438,503,313" alt="S02" href="S22.html">',
            '</map>',
        ]

    @staticmethod
    def _close_page():
        return ["</body>",
                "</html>"]

    @staticmethod
    def relative_file(filename):
        return re.sub(r"^.*/src/", "./src/", filename)

    def svg_atools_latiss(self, relative_file, page):
        """This handles image maps:"""
        page.extend(['<svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 3600 1800">',
                     f'<image width="3600" height="1800" xlink:href="{relative_file}"></image>',
                     '<a xlink:href="./c10.html"> <rect x="317" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c11.html"> <rect x="490" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c12.html"> <rect x="670" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c13.html"> <rect x="845" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c14.html"> <rect x="1020" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c15.html"> <rect x="1197" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c16.html"> <rect x="1375" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c17.html"> <rect x="1550" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c07.html"> <rect x="1550" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c06.html"> <rect x="1375" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c05.html"> <rect x="1197" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c04.html"> <rect x="1020" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c03.html"> <rect x="845" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c02.html"> <rect x="670" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c01.html"> <rect x="490" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c00.html"> <rect x="317" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '</svg>'])

    def svg_atools_comcam(self, relative_file, page):
        page.extend(['<svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 3600 1800">',
                     f'<image width="3600" height="1800" xlink:href="{relative_file}"></image>',
                     '<a xlink:href="./S00.html"> <rect x="317" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./S01.html"> <rect x="490" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./S02.html"> <rect x="670" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./S10.html"> <rect x="845" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./S11.html"> <rect x="1020" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./S12.html"> <rect x="1197" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./S20.html"> <rect x="1375" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./S21.html"> <rect x="1550" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./S22.html"> <rect x="1550" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '</svg>'])

    def svg_mosaic_latiss(self, relative_file, page):
        page.extend(['<svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 3600 1800">',
                     f'<image width="3600" height="1800" xlink:href="{relative_file}"></image>',
                     '<a xlink:href="./c10.html"> <rect x="317" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c11.html"> <rect x="490" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c12.html"> <rect x="670" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c13.html"> <rect x="845" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c14.html"> <rect x="1020" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c15.html"> <rect x="1197" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c16.html"> <rect x="1375" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c17.html"> <rect x="1550" y="216" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c07.html"> <rect x="1550" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c06.html"> <rect x="1375" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c05.html"> <rect x="1197" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c04.html"> <rect x="1020" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c03.html"> <rect x="845" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c02.html"> <rect x="670" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c01.html"> <rect x="490" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '<a xlink:href="./c00.html"> <rect x="317" y="911" fill="#fff" opacity="0" width="176" height="692"></rect> </a>',
                     '</svg>'])

    def svg_mosaic_comcam(self, relative_file, page):
        page.extend(['<svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 640 480">',
                     f'<image width="640" height="480" xlink:href="{relative_file}"></image>',
                     '<a xlink:href="./S00.html"> <rect x="111" y="306" fill="#fff" opacity="0" width="128" height="128"></rect> </a>',
                     '<a xlink:href="./S01.html"> <rect x="111" y="175" fill="#fff" opacity="0" width="128" height="128"></rect> </a>',
                     '<a xlink:href="./S02.html"> <rect x="111" y="50" fill="#fff" opacity="0" width="128" height="128"></rect> </a>',
                     '<a xlink:href="./S10.html"> <rect x="244" y="306" fill="#fff" opacity="0" width="128" height="128"></rect> </a>',
                     '<a xlink:href="./S11.html"> <rect x="244" y="175" fill="#fff" opacity="0" width="128" height="128"></rect> </a>',
                     '<a xlink:href="./S12.html"> <rect x="244" y="50" fill="#fff" opacity="0" width="128" height="128"></rect> </a>',
                     '<a xlink:href="./S20.html"> <rect x="376" y="306" fill="#fff" opacity="0" width="128" height="128"></rect> </a>',
                     '<a xlink:href="./S21.html"> <rect x="376" y="175" fill="#fff" opacity="0" width="128" height="128"></rect> </a>',
                     '<a xlink:href="./S22.html"> <rect x="376" y="50" fill="#fff" opacity="0" width="128" height="128"></rect> </a>',
                     '</svg>'])


    def image_handler(self, image_filename, page):
        """Handle images."""
        relative_file = self.relative_file(image_filename)
        do_simple = True
        if "Mosaic64" not in relative_file:
            # Skip Mosaic64
            useMap = ''
            if self.INST == 'LATISS':
                if "Mosaic" in relative_file:
                    # This is a fits_to_png thing
                    useMap = f"usemap='#mosaic_latiss'"
                else:
                    self.svg_atools_latiss(relative_file, page)
                    do_simple = False
            elif self.INST in ('LSSTComCam', 'LSSTComCamSim'):
                if "Mosaic" in relative_file:
                    # This is a fits_to_png thing
                    self.svg_mosaic_comcam(relative_file, page)
                    do_simple = False
                else:
                    useMap = f"usemap='#atools_comcam'"

            # page.append(f'<a href="{relative_file}">')
            if do_simple:
                page.append(f'<img class="center-fit" {useMap} src="{relative_file}">')
            # page.append('</a>')

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
        fig = Figure()
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
                if len(line) > 0 and line[0] == '#':
                    print("<pre class='pre-comment'>", line, "</pre>", file=ff)
                else:
                    print("<pre>", line, "</pre>", file=ff)
            for line in self._close_page():
                print(line, file=ff)


def main():
    reporter = CpvReporter()
    reporter.generate_report()
