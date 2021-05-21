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

    log.debug("CMGG - Updating config")
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
    # Overwrite default search pattern and set 'shared' to false to avoid running the module from core mqc
    ## Sampletracking
    if "sampletracking/crosscheckfingerprints" not in config.sp:
        # add new key
        config.update_dict(
            config.sp,
            {"sampletracking/crosscheckfingerprints": {"contents": "CrosscheckFingerprints", "shared": False}},
        )
        # overwrite old key
        config.update_dict(
            config.sp, {"picard/crosscheckfingerprints": {"fn": "nonexistent", "shared": False}},
        )

    ## Demultiplex
    if "picard_demultiplex/libraryparams" not in config.sp:
        config.update_dict(
            config.sp,
            {"picard_demultiplex/libraryparams": {"contents": "CollectIlluminaBasecallingMetrics", "shared": True,}},
        )

    if "picard_demultiplex/collectilluminabasecallingmetrics" not in config.sp:
        config.update_dict(
            config.sp,
            {
                "picard_demultiplex/collectilluminabasecallingmetrics": {
                    "contents": "CollectIlluminaBasecallingMetrics",
                    "shared": False,
                }
            },
        )
        config.update_dict(
            config.sp, {"picard/collectilluminabasecallingmetrics": {"fn": "nonexistent", "shared": False}},
        )

    if "picard_demultiplex/collectilluminalanemetrics" not in config.sp:
        config.update_dict(
            config.sp,
            {
                "picard_demultiplex/collectilluminalanemetrics": {
                    "contents": "CollectIlluminaLaneMetrics",
                    "shared": False,
                }
            },
        )
        config.update_dict(
            config.sp, {"picard/collectilluminalanemetrics": {"fn": "nonexistent", "shared": False}},
        )

    if "picard_demultiplex/extractilluminabarcodes" not in config.sp:
        config.update_dict(
            config.sp,
            {"picard_demultiplex/extractilluminabarcodes": {"contents": "ExtractIlluminaBarcodes", "shared": False}},
        )
        config.update_dict(
            config.sp, {"picard/extractilluminabarcodes": {"fn": "nonexistent", "shared": False}},
        )

    if "picard_demultiplex/markilluminaadapters" not in config.sp:
        config.update_dict(
            config.sp,
            {"picard_demultiplex/markilluminaadapters": {"contents": "MarkIlluminaAdapters", "shared": False}},
        )
        config.update_dict(
            config.sp, {"picard/markilluminaadapters": {"fn": "nonexistent", "shared": False}},
        )
