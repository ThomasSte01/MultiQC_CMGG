"""
The project is configured in pyproject.toml. This script is left for the editable installation.
"""

from setuptools import setup,find_packages  # type: ignore

setup(
    name='multiqc_cmgg',
    packages=find_packages(),
    include_package_data=True,
)