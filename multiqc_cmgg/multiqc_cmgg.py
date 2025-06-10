#!/usr/bin/env python
"""MultiQC example plugin functions

We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""

from __future__ import print_function

from multiqc import config
import importlib_metadata
import logging

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")


# Add default config options for the things that are used in MultiQC_NGI
def plugin_execution_start():
    """Code to execute after the config files and
    command line flags have been parsedself.

    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """

    # Plugin's version number defined in pyproject.toml:
    version = importlib_metadata.version("multiqc_cmgg")
    log.info(f"Running Example MultiQC Plugin v{version}")

    # Add to the main MultiQC config object.
    # User config files have already been loaded at this point
    #   so we check whether the value is already set. This is to avoid
    #   clobbering values that have been customised by users.

    if config.kwargs.get("disable_plugin", False):
        log.info("Plugin is disabled")
        return None

    log.debug("CMGG - Updating config")
    # Add module to module order
    config.module_order.extend(["sample_gender","coverage","MSH2_hotspot_varcount","msi_sensor_pro"])
     # Move module to the top
    config.top_modules.extend(["sample_gender","coverage","MSH2_hotspot_varcount","msi_sensor_pro"])

    list_sps=[{"sample_gender/xy": {"fn": "*_xy.tsv", "shared": False}},
            {"sample_gender/hetx": {"fn": "*_hetx.tsv", "shared": False}},
            {"sample_gender/sry": {"fn": "*_sry.tsv", "shared": False}},
            {"coverage/region_dist": {"fn": "*.mosdepth.region.dist.txt","shared": False}},
            {"MSH2_hotspot_varcount/counts": {"fn": "*.counts.txt","shared": False}},
            {"msi_sensor_pro/summary": {"fn": "*_summary_msi*","shared": False}},
            {"msi_sensor_pro/all": {"fn": "*_all_msi*","shared": False}},
    ]
    for searchp in list_sps:
        config.update_dict(config.sp,searchp)
