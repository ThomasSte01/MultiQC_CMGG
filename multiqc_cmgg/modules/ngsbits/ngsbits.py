import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc_cmgg.modules.ngsbits.sample_gender import parse_reports as samplegender_parse_reports

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The ngs-bits module parses XML output generated for several tools in the ngs-bits collection:
    * [ReadQC](https://github.com/imgag/ngs-bits/blob/master/doc/tools/ReadQC.md) for statistics on FASTQ files,
    * [MappingQC](https://github.com/imgag/ngs-bits/blob/master/doc/tools/MappingQC.md) for statistics on BAM files,
    * [SampleGender](https://github.com/imgag/ngs-bits/blob/master/doc/tools/SampleGender.md) for gender prediction based on sequencing data.
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="ngs",
            anchor="ngsbits",
            href="https://github.com/imgag/ngs-bits",
            info="Calculating statistics from FASTQ, BAM, and VCF",
            doi="10.1093/bioinformatics/btx032",
        )

        # Call submodule functions
        n = dict()
        n["sample_gender"] = samplegender_parse_reports(self)
        if n["sample_gender"] > 0:
            log.info(f"Found {n['sample_gender']} SampleGender reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
