#!/usr/bin/env python

import logging

from multiqc.utils import config
from pkg_resources import get_distribution

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.multiqc_cmgg_version = get_distribution("multiqc_cmgg").version
log.info("Running MultiQC CMGG Plugin v{}".format(config.multiqc_cmgg_version))


def update_config():
    log.info("CMGG: Updating config")
    # Nothing to do here (yet)
