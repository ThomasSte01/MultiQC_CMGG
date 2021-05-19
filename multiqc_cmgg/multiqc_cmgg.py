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
    config.module_order.append({"sampletracking": {"module_tag": ["DNA", "RNA"]}})

    # Move module to the top
    config.top_modules.append("sampletracking")

    # Disable module to avoid duplicate data
    disabled_modules = []
    for module in disabled_modules:
        del config.avail_modules[module]

    # Update search patterns
    # Overwrite default crosscheckfingerprints search pattern and set 'shared' to false to avoid running the module from core mqc
    if "sampletracking/crosscheckfingerprints" not in config.sp:
        config.update_dict(
            config.sp,
            {"sampletracking/crosscheckfingerprints": {"contents": "CrosscheckFingerprints", "shared": False}},
        )
