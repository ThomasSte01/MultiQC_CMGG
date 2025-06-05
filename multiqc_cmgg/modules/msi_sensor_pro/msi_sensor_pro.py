import logging
from collections import defaultdict
from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.utils.util_functions import update_dict
from multiqc.plots import table,bargraph
from typing import Dict, Union, List, Optional
from collections import OrderedDict

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
             name="msisensor-pro",
             info="This table show a summary of the output of the msisensor-pro pro command for MSI detection.",
        )

        # Load configurable thresholds (default to 300 if not specified)
        self.coverage_threshold = getattr(config, "msi_sensor_pro_coverage_threshold")
        self.min_sites_threshold = getattr(config, "msi_sensor_pro_min_sites")
        self.msi_high_threshold = getattr(config, "msi_high_threshold")

        # Parsing and loading data from msiSensorPro summary and all files
        data_dicts_summary = self.parse_summary()
        data_dicts_all = self.parse_all()
        msisensorpro_data, all_zero = self.prepare_msisensorpro_data(data_dicts_summary)

        # Table configuration
        config_table ={
            "id": "msi_summary",
            "title": "msi_summary",
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

        # headers for all table
        all_loci = set()
        for sample in data_dicts_all.values():
            all_loci.update(sample.keys())

        headers2 = {}
        for locus in sorted(all_loci):
            headers2[locus] = {
                "title": locus,
                "description": f"MSI status at {locus}"
            }

        # summary table
        self.add_section(
            plot=table.plot(data=data_dicts_summary, headers=headers, pconfig=config_table),
        )

        # all table
        self.add_section(
            name="msisensor-pro - All Loci",
            anchor="msisensor-pro - all loci",
            description="Detailed MSI status for every loci with the coverage in brackets.",
            plot=table.plot(
                data=data_dicts_all,
                headers=headers2,
                pconfig={
                    "id": "msiSensorPro_all_table",
                    "title": "msiSensorPro - All Site Metrics"
                }
            )
        )
        if not all_zero:
            categories = {
                "MSS": {
                    "name": "MSS",
                    "color": "#2ecc71"
                },
                "MSI-high": {
                    "name": "MSI-high", 
                    "color": "#e74c3c"
                },
                "Low-coverage": {
                    "name": "Low-coverage",
                    "color": "#f39c12"
                }
            }

            # Bargraph configuration
            self.add_section(
                name="msisensor-pro - Bargraph",
                anchor="msisensor-pro - Bargraph",
                description="This graph visualizes the MSI status per sample in a bargraph.",
                plot=bargraph.plot(
                    data=msisensorpro_data,
                    cats=categories,
                    pconfig={
                        "id": "msiSensorPro_bargraph",
                        "title": "MSI Sensor Pro Summary",
                        "ylab": "Percentage of unstable sites",
                        "ymin": 0,
                        "ymax": 100,
                        "cpswitch_counts_label": "Number of Sites",
                    }
                )
            )
        else:
            log.info("Skipping bargraph: All samples have 0% unstable sites")

    def prepare_msisensorpro_data(self, data_summary):
        """
        Transform summary data to msisensorpro score with MSI classification
        """
        msisensorpro_data = {}
        all_zero = True
        
        for sample_name, sample_data in data_summary.items():

            msisensorpro_score = sample_data["perc"]
            
            # Classify MSI status based on thresholds
            if sample_data["num_sites"] <= self.min_sites_threshold:
                msi_status = "Low-coverage"
            elif msisensorpro_score >= self.msi_high_threshold:  # 30% threshold for MSI-high
                msi_status = "MSI-high"
            else:
                msi_status = "MSS"
            
            # Structure data for bargraph - each sample gets assigned to one category
            sample_entry = {
                "MSS": msisensorpro_score if msi_status == "MSS" else 0,
                "MSI-high": msisensorpro_score if msi_status == "MSI-high" else 0,
                "Low-coverage": msisensorpro_score if msi_status == "Low-coverage" else 0
            }
            msisensorpro_data[sample_name] = sample_entry
            if any(value != 0.0 for value in sample_entry.values()):
                all_zero = False
        
        return msisensorpro_data, all_zero
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
    
    def parse_all(self) -> Dict[str, Dict]:
        """
        Parse the msiSensorPro all file into a loci-centric structure.
        """
        sample_data = defaultdict(dict)
        all_loci = set()

        for f in self.find_log_files("msi_sensor_pro/all", filecontents=True, filehandles=False):
            s_name = self.clean_s_name(f["fn"], f)
            lines = f["f"].splitlines()
            
            for line in lines[1:]:  # Skip header
                parts = line.strip().split("\t")
                if len(parts) < 10:
                    log.warning(f"Skipping line in {s_name} due to insufficient columns")
                    continue

                chrom, pos = parts[0], parts[1]
                locus_id = f"{chrom}:{pos}"
                all_loci.add(locus_id)
                pro_p = float(parts[6])
                coverage = int(parts[8])
                threshold = float(parts[9])

                # Classify MSI status
                if pro_p > threshold:
                    status = f"Unstable ({coverage})"
                elif coverage < self.coverage_threshold:
                    status = f"Low-coverage ({coverage})"
                else:
                    status = f"Stable ({coverage})"

                sample_data[s_name][locus_id] = status
        
        for sample in sample_data.values():
            missing = all_loci - set(sample.keys())
            for locus in missing:
                sample[locus] = "Low-coverage (0)"
        
        return sample_data

