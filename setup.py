#!/usr/bin/env python
"""
MultiQC_CMGG is a plugin for MultiQC containing customized modules and templates.
"""

from setuptools import setup, find_packages

version = "0.0.1"

setup(
    name="multiqc_cmgg",
    version=version,
    author="Matthias De Smet",
    author_email="11850640+matthdsm@users.noreply.github.com",
    description="MultiQC plugin for interal use at CMGG",
    long_description=__doc__,
    keywords="bioinformatics",
    url="https://github.com/CenterForMedicalGeneticsGhent/MultiQC_CMGG",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=["multiqc>=1.10"],
    entry_points={
        "multiqc.hooks.v1": ["config_loaded = multiqc_cmgg.multiqc_cmgg:update_config",],
        "multiqc.modules.v1": [
            "sampletracking = multiqc_cmgg.modules.sampletracking.sampletracking:MultiqcModule",
            "picard_demultiplex = multiqc_cmgg.modules.picard_demultiplex.demultiplex:MultiqcModule",
        ],
        "multiqc.templates.v1": ["cmgg = multiqc_cmgg.templates.cmgg",],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Environment :: Web Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)

