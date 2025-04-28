#!/usr/bin/env python

import logging

from multiqc.utils import config
from pkg_resources import get_distribution

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.multiqc_cmgg_version = get_distribution("multiqc_cmgg").version
log.info(f"Running MultiQC CMGG Plugin v{config.multiqc_cmgg_version}")


def update_config() -> None:
    """
    Update MultiQC config object
    * Update module order
    * Disable unnecessary modules to avoid duplicate data
    * Update search patterns
    """

    if config.kwargs.get("disable_plugin", True):
        log.info("disable_plugin was True")
        return None

    log.debug("CMGG - Updating config")
    # Add module to module order
    config.module_order.extend(
        [
            {"ngsbits": {"module_tag": ["DNA"]}},
            # {"coverage": {"module_tag": ["DNA", "RNA"]}},
        ]
    )

    # Move module to the top
    config.top_modules.extend(["ngsbits", "coverage"])

    # Disable module to avoid duplicate data
    disabled_modules = []
    for module in disabled_modules:
        del config.avail_modules[module]

    # Update search pattern
    update_search_patterns()

    # Update fn cleanup
    update_fn_cleanup()


def update_search_patterns() -> None:
    """
    Update search patterns
    Overwrite default search pattern and set 'shared' to false to avoid running the module from core mqc
    :return: None
    """
    ##Sample_gender
    if "ngsbits/sample_gender" not in config.sp:
        # add new key
        config.update_dict(
            config.sp,
            {"ngsbits/sample_gender": {"contents": "SampleGenders", "shared": False}},
        )
        # overwrite old key
        config.update_dict(
            config.sp, {"ngsbits/sample_gender": {"fn": "nonexistent", "shared": False}},
        )

def update_fn_cleanup() -> None:
    """
    Update filename cleanup patterns
    :return: None
    """
    config.fn_clean_exts.extend(["_duplicate_metrics", "_samtools_stats", "_samtools_idxstats", ".unaligned"])
