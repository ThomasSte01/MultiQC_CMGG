import logging
from collections import defaultdict
from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.utils.util_functions import update_dict
from multiqc.plots import table
from typing import Dict, Union, List, Optional

log = logging.getLogger(__name__)

def genstats_cov_thresholds(cum_fraction_by_cov: Dict[int, float], threshs: List[int]) -> Dict[str, float]:
    genstats: Dict[str, float] = {}
    for t in threshs:
        genstats[f"{t}_x_pc"] = cum_fraction_by_cov.get(t, 0.0) * 100.0
    return genstats

def calc_median_coverage(cum_fraction_by_cov) -> Optional[float]:
    median_cov = None
    for this_cov, cum_fraction in sorted(cum_fraction_by_cov.items(), reverse=True):
        if cum_fraction >= 0.5:
            median_cov = this_cov
            break
    return median_cov

class MultiqcModule(BaseMultiqcModule):
    """
    Mosdepth can generate several output files all with a common prefix and different endings:

    - per-base depth (`{prefix}.per-base.bed.gz`),
    - mean per-window depth given a window size (`{prefix}.regions.bed.gz`, if a BED file provided with `--by`),
    - mean per-region given a BED file of regions (`{prefix}.regions.bed.gz`, if a window size provided with `--by`),
    - a distribution of proportion of bases covered at or above a given threshhold for each chromosome and genome-wide (`{prefix}.mosdepth.global.dist.txt` and `{prefix}.mosdepth.region.dist.txt`),
    - quantized output that merges adjacent bases as long as they fall in the same coverage bins (`{prefix}.quantized.bed.gz`),
    - threshold output to indicate how many bases in each region are covered at the given thresholds (`{prefix}.thresholds.bed.gz`)
    - summary output providing region length, coverage mean, min, and max for each region. (`{prefix}.mosdepth.summary.txt`)

    The MultiQC module plots coverage distributions from 1 output:

    - `{prefix}.mosdepth.region.dist.txt`

    Plotting the percentage of the genome covered at a threshold in the General Stats section.
    The default thresholds are 20,30 which can be customised in the config as follows:

    ```yaml
    mosdepth_config:
      general_stats_coverage:
        - 20
        - 30
    ```
    ```yaml
    general_stats_coverage_hidden:
        - 30
    ```
    """
    # Initialise the parent module Class object
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="coverage",
            anchor="coverage",
            href="https://github.com/brentp/mosdepth",
            info="Fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing",
            doi="10.1093/bioinformatics/btx699",
        )     
        # Load in thresholds and hidden thresholds from config
        threshs, hidden_threshs = config.get_cov_thresholds("coverage_config")
        thresh = threshs[0]
        
        genstats_by_sample: Dict[str, Dict[str, Union[int, float]]] = defaultdict(dict)  # mean coverage

        # Filter out any samples from --ignore-samples
        data_dicts_region = self.parse_cov_dist()
    
        data_dicts_region = self.ignore_samples(data_dicts_region)

        samples_region = set(data_dicts_region.keys())
        samples_found = samples_region
        if not samples_found:
            raise ModuleNoSamplesFound
        log.info(f"Found reports for {len(samples_found)} samples")
        
        # Adding hide buttons to report, based on run names
        list_yaml_configs=["show_hide_buttons","show_hide_mode","show_hide_patterns,show_hide_color"]

        for yaml_header in list_yaml_configs:
            if not hasattr(config, yaml_header) or getattr(config,yaml_header) is None:
                setattr(config,yaml_header,{})
        list_show_hide_mode=[]
        list_show_hide_names=[]
        list_show_hide_color=[]
        more_less90={}

        for f in self.find_log_files(f"coverage/region_dist", filecontents=False, filehandles=False):
            s_name = self.clean_s_name(f["fn"], f)
            
            name=s_name.split("_",1)[0] if "_" in s_name else s_name
            if name not in list_show_hide_names:
                    list_show_hide_names.append(name)
                    list_show_hide_mode.append("show")
    
            if "_" not in s_name and s_name in data_dicts_region:
                if name not in more_less90:
                    more_less90[name]=""
                if data_dicts_region[s_name][f'{str(thresh)}_x_pc']<90:
                    more_less90[name]="#d9534f"

                if data_dicts_region[s_name][f'{str(thresh)}_x_pc']>= 90 and more_less90[name]!="#d9534f":
                    more_less90[name]="#5cb85c"
        
        config.show_hide_buttons=list_show_hide_names
        config.show_hide_patterns=list_show_hide_names
        config.show_hide_mode=list_show_hide_mode
        for color in more_less90.values():
            list_show_hide_color.append(color)
        config.show_hide_color=list_show_hide_color

        # Generalstat_headers setup
        genstats_headers = {}

        
        for t in threshs:
            genstats_headers[f"{t}_x_pc"] = {
                "title": f"≥ {t}X",
                "description": f"Fraction of genome with at least {t}X coverage",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": False,
                "hidden": t in hidden_threshs,
                "cond_formatting_rules":{"pass":[{"gt": 90},{"eq": 90}],"fail":[{"lt":90}]},
                "cond_formatting_colours":[{"pass":"#5cb85c"},{"fail":"#d9534f"}],
            }
        # Add coverage summary to General Stats
        genstats_headers.update(
            {
                "median_coverage": {
                    "title": "Median",
                    "description": "Median coverage",
                    "min": 0,
                    "suffix": "X",
                    "scale": "BuPu",
                },
            },
        )
        self.general_stats_addcols(data_dicts_region, genstats_headers)

    def parse_cov_dist(self,) -> Dict[str, Dict[str, Union[float, int, None]]]:
        # per chromosome average coverage
        genstats_by_sample: Dict[str, Dict[str, Union[float, int, None]]] = dict()
        threshs, hidden_threshs = config.get_cov_thresholds("coverage_config")
        for f in self.find_log_files(f"coverage/region_dist", filecontents=False, filehandles=True):
            s_name = self.clean_s_name(f["fn"], f)
            self.add_data_source(f, s_name=s_name, section="genome_results")
            cum_fraction_by_cov: Dict[int, float] = dict()
            for line in f["f"]:
                contig, cutoff_reads, bases_fraction = str(line).split("\t")
                if bases_fraction == "0.00\n":
                    continue
                # Parse cumulative coverage
                if contig == "total":
                    cum_fraction_by_cov[int(cutoff_reads)] = float(bases_fraction)
            genstats_by_sample[s_name] = {}
            for k, v in genstats_cov_thresholds(cum_fraction_by_cov, threshs).items():
                genstats_by_sample[s_name][k] = v
            genstats_by_sample[s_name]["median_coverage"] = calc_median_coverage(cum_fraction_by_cov)
        
        return genstats_by_sample
