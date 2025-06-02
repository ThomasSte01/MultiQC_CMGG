import logging
from collections import defaultdict
from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.utils.util_functions import update_dict
from multiqc.plots import table,bargraph
from typing import Dict, Union, List, Optional

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
             name="msi_sensor_pro",
             info="This table show the results of msiSensorPro for msi detection and some metrics.",
        )
        # Parsing and loading data from msiSensorPro summary and all files
        data_dicts_summary = self.parse_summary()
        data_dicts_all = self.parse_all()

        # Table configuration
        config_table ={
            "id": "table_name",
            "title": "title",
        }
        headers={
            "num_sites": {
                "title": "Number of sites",
            },
            "num_unstable_sites": {
                "title": "Number of unstable sites",
            },
            "perc": {
                "title": "Percentage of unstable sites",
                "format": "{:.2f}",
                "suffix": "%",
            },

        }

        self.add_section(
            plot=table.plot(data=data_dicts_summary, headers=headers, pconfig=config_table),
        )
    
        # Bargraph configuration
        self.add_section(
            name="msiSensorPro",
            anchor="msiSensorPro",
            description="This section contains the results of msiSensorPro for msi detection and some metrics.",
            plot=bargraph.plot(
                data=data_dicts_summary,
                pconfig={
                    "id": "msiSensorPro_bargraph",
                    "title": "MSI Sensor Pro Summary",
                    "ylab": "Percentage of unstable sites",
                    "ymin": 0,
                    "ymax": 100,
                }
            )
        )
    # Parsing summary file for msiSensorPro
    def parse_summary(self,):
        """
        Parse the msiSensorPro summary file.
        """
        data_summary: Dict[str, Dict[str, float]] = {}
        for f in self.find_log_files("msi_sensor_pro/summary", filecontents=True, filehandles=False):
            s_name = self.clean_s_name(f["fn"], f)
            lines=f["f"].splitlines()
            header = lines[0]
            for line in lines:
                if line != header:
                    num_sites, num_unstable_sites, perc = line.split("\t")
                    data_summary[s_name] = {
                    "num_sites": int(num_sites),
                    "num_unstable_sites": int(num_unstable_sites),
                    "perc": float(perc)
                    }
            log.info(data_summary)
        return data_summary
    
    # Parsing all file for msiSensorPro
    def parse_all(self,) -> Dict[str, Dict[str, List]]:
        """
        Parse the msiSensorPro all file.
        """
        data_all: Dict[str, Dict[str, List[str,int,float ]]] = defaultdict(dict)
        for f in self.find_log_files("msi_sensor_pro/all", filecontents=True, filehandles=False):
            s_name = self.clean_s_name(f["fn"], f)
            lines=f["f"].splitlines()
            header = lines[0]
            # log.info(lines)
            for line in lines:
                parts = str(line).split("\t")
                if len(parts) < 10:
                    log.warning(f"Skipping line in {s_name} due to insufficient  filled in columns")
                    continue
                if line != header:
                # Might have to change this if it causes issues for handeling of table
                    data_all[s_name][f"{parts[0]}_{parts[1]}"] = {
                        "chrom": parts[0],
                        "loc": int(parts[1]),
                        "left_flank_bases": parts[2],
                        "repeat_times": int(parts[3]),
                        "repeat_unit_bases": parts[4],
                        "right_flank_bases": parts[5],
                        "pro_p": float(parts[6]),
                        "pro_q": float(parts[7]),
                        "CovReads": int(parts[8]),
                        "threshold": float(parts[9]),
                    }
        
        return data_all
        
