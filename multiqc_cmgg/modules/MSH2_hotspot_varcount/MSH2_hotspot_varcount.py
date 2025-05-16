import logging
from multiqc import config
from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table
from typing import Dict, Union

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")

log.info("started MSH2_hotspot_varcount")
log.info(f"threshold Sanger has been set to {config.MSH2_hotspot_varcount_config["sanger_threshold"]}")
class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        
        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
          name="MSH2_hotspot_varcount",
          info="Analysis module used for MSH2 variance counting",
        )

        # Find and load any input files for this module
        MSH2_varcount_data : Dict[str, Dict[str, Union[float, str]]] = dict()
        
        for f in self.find_log_files("MSH2_hotspot_varcount/counts"):
            self.add_data_source(f)
            s_name = f["s_name"]
            parsed =parse_file(f["f"])
            if s_name not in  MSH2_varcount_data:
                MSH2_varcount_data[s_name] = parsed
            
            if s_name in  MSH2_varcount_data:
                MSH2_varcount_data[s_name].update(parsed)
            # Filter to strip out ignored sample names
                MSH2_varcount_data = self.ignore_samples( MSH2_varcount_data)

        # Debug for amount of reports found
        n_reports_found = len(MSH2_varcount_data)
        if n_reports_found > 0:
            log.debug(f"Found {len(MSH2_varcount_data)} MSH2_hotspot_varcount reports")

        if n_reports_found == 0:
            log.debug("No MSH2_hotspot_varcount reports found")

        # Write parsed report data to a file
        self.write_data_file(MSH2_varcount_data, "multiqc_MSH2_hotspot_varcount")
        self.add_software_version(None)

        # log.info(MSH2_varcount_data)
        # Add MSH2_hotspot Table
        config_table = {
            "id": "MSH2_hotspot_varcount",
            "title":"MSH2_hotspot_varcount",
            "sort_rows":True,
            "no_violin": True,
        }

        headers = {
            "MSH2_c.942+3_wt":{
                "title" :"WT count",
                "description":"counts of WT",
                "scale":"PuBu",
            },
            "MSH2_c.942+3A>T":{
                "title" : "MSH2_c.942+3A>T (counts)",
                "description":"frequency of 3A>T and counts",
                "cond_formatting_rules":{"sanger":[{"s_contains":"- Sangeren!"}]},
                "cond_formatting_colours":[{"sanger":"#EE4B2B"}],
                "scale":False,
            },
            "MSH2_c.942+2T>A":{
                "title" : "MSH2_c.942+2T>A (counts)",
                "description":"frequency of 2T>A and counts",
                "cond_formatting_rules":{"sanger":[{"s_contains":"- Sangeren!"}]},
                "cond_formatting_colours":[{"sanger":"#EE4B2B"}],
                "scale":False,
            },
            "MSH2_c.942+2T>C":{
                "title" : "MSH2_c.942+2T>C (counts)",
                "description":"frequency of 2T>C and counts",
                "cond_formatting_rules":{"sanger":[{"s_contains":"- Sangeren!"}]},
                "cond_formatting_colours":[{"sanger":"#EE4B2B"}],
                "scale":False,
            },
            "MSH2_c.942+2T>G":{
                "title" : "MSH2_c.942+2T>G (counts)",
                "description":"frequency of 2T>G and counts",
                "cond_formatting_rules":{"sanger":[{"s_contains":"- Sangeren!"}]},
                "cond_formatting_colours":[{"sanger":"#EE4B2B"}],
                "scale":False,
            },
            
        }
        self.add_section(
            plot=table.plot(data=MSH2_varcount_data, headers=headers, pconfig=config_table),
        )
def parse_file(f: str) -> Dict[str, Union[float, str]]:
    """
    Parses a single samplegender TSV file content and returns a dictionary
    with the relevant data from columns 2-6.
    """
    parsed_data: Dict[str, Union[float, str]] = {}
    lines = f.splitlines()
    
    if len(lines) < 3:
        # Not enough data, return an empty dictionary
        return parsed_data
    headers = lines[2].strip().split(" ")[1:6]
    values = lines[3].strip().split(" ")[1:6]
 
    for key, value in zip(headers, values):
        parsed_data[key] = value

    #calculating frequency of mutation:
    for variant,counts in parsed_data.items():
        if variant != "MSH2_c.942+3_wt":
            freq=round((int(counts)/(int(parsed_data["MSH2_c.942+3_wt"])+int(counts)))*100 ,2)
            # Determining need for sanger
            if freq >= config.MSH2_hotspot_varcount_config["sanger_threshold"]:
                parsed_data[variant]=f"{freq}({counts}) - Sangeren!"
            else:
                parsed_data[variant]=f"{freq}({counts})"
        else:
            parsed_data[variant]=counts

    return parsed_data