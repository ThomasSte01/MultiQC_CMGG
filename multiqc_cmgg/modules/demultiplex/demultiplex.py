#!/usr/bin/env python

""" MultiQC module to parse output from Picard demultiplexing tools"""

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc.modules.base_module import BaseMultiqcModule

# Import the Demux submodules
from . import CollectIlluminaBasecallingMetrics
from . import CollectIlluminaLaneMetrics
from . import ExtractIlluminaBarcodes
from . import MarkIlluminaAdapters

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """Picard is a collection of scripts. This MultiQC module
    supports some but not all. The code for each script is split
    into its own file and adds a section to the module output if
    logs are found."""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Picard Demultiplexing", anchor="picard", info="Subset of Picard tools used for demultiplexing",
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()

        # Call submodule functions
        n["CollectIlluminaBasecallingMetrics"] = CollectIlluminaBasecallingMetrics.parse_reports(self)
        if n["CollectIlluminaBasecallingMetrics"] > 0:
            log.info(
                "Found {} CollectIlluminaBasecallingMetrics reports".format(n["CollectIlluminaBasecallingMetrics"])
            )

        n["CollectIlluminaLaneMetrics"] = CollectIlluminaLaneMetrics.parse_reports(self)
        if n["CollectIlluminaLaneMetrics"] > 0:
            log.info("Found {} CollectIlluminaLaneMetrics reports".format(n["CollectIlluminaLaneMetrics"]))

        n["ExtractIlluminaBarcodes"] = ExtractIlluminaBarcodes.parse_reports(self)
        if n["ExtractIlluminaBarcodes"] > 0:
            log.info("Found {} ExtractIlluminaBarcodes reports".format(n["ExtractIlluminaBarcodes"]))

        n["MarkIlluminaAdapters"] = MarkIlluminaAdapters.parse_reports(self)
        if n["MarkIlluminaAdapters"] > 0:
            log.info("Found {} MarkIlluminaAdapters reports".format(n["MarkIlluminaAdapters"]))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise UserWarning

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)


def read_histogram(self, program_key, program_name, headers, formats):
    """
    Reads a Picard HISTOGRAM file.

    Args:
        self: the Picard QC module
        program_key: the key used to find the program (ex. picard/quality_by_cycle)
        program_name: the program key in the header to find the I/INPUT line
        headers: the list of expected headers for the histogram
        formats: the list of methods to apply to re-format each field (on a given row)
    """
    all_data = OrderedDict()

    assert len(formats) == len(headers)

    # Go through logs and find Metrics
    for f in self.find_log_files(program_key, filehandles=True):
        lines = iter(f["f"])

        # read through the header of the file to obtain the
        # sample name
        clean_fn = lambda n: self.clean_s_name(n, f["root"])
        s_name = read_sample_name(lines, clean_fn, program_name)
        if s_name is None:
            continue

        sample_data = OrderedDict()

        try:
            # skip to the histogram
            line = next(lines)
            while not line.startswith("## HISTOGRAM"):
                line = next(lines)

            # check the header
            line = next(lines)
            if headers != line.strip().split("\t"):
                continue

            # slurp the data
            line = next(lines).rstrip()
            while line:
                fields = line.split("\t")
                assert len(fields) == len(headers)
                for i in range(len(fields)):
                    fields[i] = formats[i](fields[i])

                sample_data[fields[0]] = OrderedDict(zip(headers, fields))
                line = next(lines).rstrip()

        except StopIteration:
            pass

        # append the data
        if sample_data:
            all_data[s_name] = sample_data

    return self.ignore_samples(all_data)
