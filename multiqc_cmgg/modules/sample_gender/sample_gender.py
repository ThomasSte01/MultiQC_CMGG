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
             name="Sex prediction",
             info="This table show the results of the sex prediction with the expected sex and some metrics.",
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
        # log.info(f" dict with results: {samplegender_data}")
        
        # Calculating certainty of gender based of amount of times gender was determined and calculated sex
        list_gender_methods=['gender_xy','gender_hetx','gender_sry']
        for sample in samplegender_data.keys():
            count_M=0
            count_F=0
            for method in list_gender_methods:
                if samplegender_data[sample][method]=="M":
                    count_M+=1 
                if samplegender_data[sample][method]=="F":
                    count_F+=1
            if count_M >= count_F:
                percentage=abs(count_M/3)
                calc_sex="M"
            elif count_F >= count_M:
                percentage=abs(count_F/3)
                calc_sex="F"    
            samplegender_data[sample].update({"certainty":percentage})
            samplegender_data[sample].update({"calc_gender":calc_sex})

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
            "id": "Sex prediction",
            "title":"Sex prediction",
        }

        headers = {
            "certainty": {
                "title": "Certainty of sex match",
                "description":"match_certainty",
                "format": "{:.0%}",
                "scale":False,
                "max":1,
                "cond_formatting_rules":{"pass":[{"eq": 1}],"warn":[{"lt":1}],"fail":[{"lt":0.4}]},
                "cond_formatting_colours":[{"pass":"#5cb85c"},{"warn":"#f0ad4e"}, {"fail":"#d9534f"}],
            },
            "calc_gender":{
                "title": "Calculated Sex",
                "description":"calculated_sex",
                "scale":False,
            },
            "gender_xy": {
                "title": "Calculated sex (XY method)",
                "description": "The predicted gender based on chromosome read ratios.",
                "scale": False,
            },
            "gender_sry": {
                "title": "Calculated sex (SRY method)",
                "description": "The predicted gender coverage of SRY gene.",
                "scale": False,
            },            
            "gender_hetx": {
                "title": "Calculated sex (HETX method)",
                "description": "The predicted gender fraction of heterozygous variants on X chromosome.",
                "scale": False,
                "cond_formatting_rules":{"unknown":[{"s_eq":"unknown (too few SNPs)"}]},
                "cond_formatting_colours":[{"unknown":"#808080"}],
            },
            "ratio_chry_chrx": {
                "title": "ChrY/ChrX reads ratio",
                "description": "The ratio of reads mapped to ChrY vs ChrX.",
                "min": 0,
                "format": "{:.4f}",
                "scale": "Purples",
            },
            "coverage_sry": {
                "title": "Coverage SRY",
                "description": "Coverage of SRY in chrY (SRY)",
                "min": 0,
                "scale": "Blues",
                "format": "{:,.2f}",
            },
            "het_fraction": {
                "title": "Fraction HETX",
                "description": "Fraction of heterozygous SNPs in chrX (HETX)",
                "min": 0,
                "scale": "Reds",
                "format": "{:,.4f}",
            },
        }

        self.add_section(
            plot=table.plot(data=samplegender_data, headers=headers, pconfig=config_table),
        )
    
def parse_file(f: str) -> Dict[str, Union[float, str]]:
    """
    Parses a single samplegender TSV file content and returns a dictionary
    with the relevant data from columns 2-6.
    """
    parsed_data: Dict[str, Union[float, str]] = {}
    lines = f.splitlines()
    
    if len(lines) < 2:
        # Not enough data, return an empty dictionary
        return parsed_data
    headers = lines[0].strip().split("\t")[1:6]
    values = lines[1].strip().split("\t")[1:6]

    #Changing gender names to abbriviations
    if values[0]=="male":
        values[0]="M"
    elif values[0]=="female":
        values[0]="F"

    # Changes gender so that every method result can be separated
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