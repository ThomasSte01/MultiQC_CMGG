import logging
from multiqc import config
from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table
from typing import Dict, Union
import math


# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

log.info("started module")
class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name="sample_gender",
            anchor="ngsbits_sample_gender",
            href="https://github.com/imgag/ngs-bits",
            info="Calculating statistics from FASTQ, BAM, and VCF",
            doi="10.1093/bioinformatics/btx032",
        )

        # Find and load any input files for this module
        samplegender_data : Dict[str, Dict[str, Union[float, str]]] = dict()
        # log.info(f"find log files: {list(self.find_log_files('sample_gender/xy'))}")
        method_dict={"sample_gender/xy":"_xy","sample_gender/hetx":"_hetx","sample_gender/sry":"_sry"}
        for method,extension in method_dict.items():

            for f in self.find_log_files(method):
                self.add_data_source(f)
                s_name = f["s_name"]
                s_name = s_name.replace(extension,"")
                parsed =parse_file(f["f"])

                if s_name not in samplegender_data:
                    samplegender_data[s_name] = parsed
                
                if s_name in samplegender_data:
                    samplegender_data[s_name].update(parsed)
                # Filter to strip out ignored sample names
                samplegender_data = self.ignore_samples(samplegender_data)
        log.debug(f" dict with results: {samplegender_data}")
        
        n_reports_found = len(samplegender_data)
        if n_reports_found > 0:
            log.debug(f"Found {len(samplegender_data)} SampleGender reports")

        if n_reports_found == 0:
            log.debug("No SampleGender reports found")
        
        # Write parsed report data to a file
        self.write_data_file(samplegender_data, "multiqc_ngsbits_samplegender")
        self.add_software_version(None)

        # Add samplegender Table
        config_table = {
            "id": "samplegender",
            "title": "samplegender",
        }

        headers = {
            "gender_xy": {
                "title": "Calculated sex (XY method)",
                "description": "The predicted gender based on chromosome read ratios.",
                "namespace": "ngsbits",
                "scale": False,
            },
            "gender_sry": {
                "title": "Calculated sex (SRY method)",
                "description": "The predicted gender coverage of SRY gene.",
                "namespace": "ngsbits",
                "scale": False,
            },            
            "gender_hetx": {
                "title": "Calculated sex (HETX method)",
                "description": "The predicted gender fraction of heterozygous variants on X chromosome.",
                "namespace": "ngsbits",
                "scale": False,
            },
            "ratio_chry_chrx": {
                "title": "ChrY/ChrX reads ratio",
                "description": "The ratio of reads mapped to ChrY vs ChrX.",
                "namespace": "ngsbits",
                "min": 0,
                "format": "{:.4f}",
                "scale": "Purples",
            },
            "coverage_sry": {
                "title": "Coverage SRY",
                "description": "Coverage of SRY in chrY (SRY)",
                "namespace": "ngsbits",
                "min": 0,
                "scale": "Blues",
                "format": "{:,.2f}",
            },
            "het_fraction": {
                "title": "Fraction HETX",
                "description": "Fraction of heterozygous SNPs in chrX (HETX)",
                "namespace": "ngsbits",
                "min": 0,
                "scale": "Reds",
                "format": "{:,.4f}",
            },
        }

        self.add_section(
            name="samplegender",
            anchor="ngsbits-samplegender",
            description='<a href="https://github.com/imgag/ngs-bits/blob/master/doc/tools/SampleGender.md" target="_blank">SampleGender</a>'
            " determines the gender of a sample from the BAM/CRAM file.",
            plot=table.plot(data=samplegender_data, headers=headers, pconfig=config_table),
        )

        for header in headers.values():
            header["hidden"] = True
            headers["gender_xy"]["hidden"] = False
            headers["gender_hetx"]["hidden"] = False
            headers["gender_sry"]["hidden"] = False

        self.general_stats_addcols(samplegender_data, headers)
    
def parse_file(f: str) -> Dict[str, Union[float, str]]:
    """
    Parses a single samplegender TSV file content and returns a dictionary
    with the relevant data from columns 2-5.
    """
    parsed_data: Dict[str, Union[float, str]] = {}
    lines = f.splitlines()
    # log.info(lines)
    if len(lines) < 2:
        # Not enough data, return an empty dictionary
        return parsed_data
    headers = lines[0].strip().split("\t")[1:6]
    values = lines[1].strip().split("\t")[1:6]

    #changes gender so that every method result can be seperated
    paramdict={"reads_chry":"gender_xy","het_fraction":"gender_hetx","coverage_sry":"gender_sry"}
    for param,test in paramdict.items():
        if param in headers:
            headers[0]=test

    for key, value in zip(headers, values):
        try:
            parsed_data[key] = float(value)
            if math.isnan(float(value)):
                parsed_data[key] = "N/A"
        except ValueError:
            parsed_data[key] = value
    return parsed_data