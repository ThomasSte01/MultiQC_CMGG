#!/usr/bin/env python

import logging

from multiqc.utils import config
from pkg_resources import get_distribution

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.multiqc_cmgg_version = get_distribution("multiqc_cmgg").version
log.info("Running MultiQC CMGG Plugin v{}".format(config.multiqc_cmgg_version))


def update_config() -> None:
    """
    Update MultiQC config object
    * Update module order
    * Disable unnecessary modules to avoid duplicate data
    * Update search patterns
    """

    log.debug("Updating config")
    # Add module to module order
    config.module_order.extend(
        [
            {"sampletracking": {"module_tag": ["DNA", "RNA"]}},
            {"demultiplex": {"module_tag": ["DNA", "RNA", "Demultiplex"]}},
        ]
    )

    # Move module to the top
    config.top_modules.extend(["sampletracking", "demultiplex"])

    # Disable module to avoid duplicate data
    disabled_modules = []
    for module in disabled_modules:
        del config.avail_modules[module]

    # Update search patterns
    ## Sampletracking
    # Overwrite default crosscheckfingerprints search pattern and set 'shared' to false to avoid running the module from core mqc
    if "sampletracking/crosscheckfingerprints" not in config.sp:
        config.update_dict(
            config.sp,
            {"sampletracking/crosscheckfingerprints": {"contents": "CrosscheckFingerprints", "shared": False}},
        )
    ## Demultiplex
    if "demux/collectilluminabasecallingmetrics" not in config.sp:
        config.update_dict(
            config.sp,
            {
                "demux/collectilluminabasecallingmetrics": {
                    "contents": "CollectIlluminaBasecallingMetrics",
                    "shared": False,
                }
            },
        )
    if "demux/collectilluminalanemetrics" not in config.sp:
        config.update_dict(
            config.sp,
            {"picard/collectilluminalanemetrics": {"contents": "CollectIlluminaLaneMetrics", "shared": False}},
        )
    if "demux/extractilluminabarcodes" not in config.sp:
        config.update_dict(
            config.sp, {"picard/extractilluminabarcodes": {"contents": "ExtractIlluminaBarcodes", "shared": False}},
        )
    if "demux/markilluminaadapters" not in config.sp:
        config.update_dict(
            config.sp, {"picard/markilluminaadapters": {"contents": "MarkIlluminaAdapters", "shared": False}},
        )

