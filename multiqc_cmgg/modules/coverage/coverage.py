import fnmatch
import logging
from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Union

from multiqc import Plot, config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.modules.qualimap.QM_BamQC import genome_fraction_helptext
from multiqc.plots import bargraph, linegraph
from multiqc.plots.linegraph import smooth_array
from multiqc.utils.util_functions import update_dict

log = logging.getLogger(__name__)


def read_config():
    cfg = getattr(config, "mosdepth_config", dict())
    if not isinstance(cfg, dict):
        return {}

    cfg["include_contigs"] = cfg.get("include_contigs", [])
    if not isinstance(cfg["include_contigs"], list):
        cfg["include_contigs"] = []

    cfg["exclude_contigs"] = cfg.get("exclude_contigs", [])
    if not isinstance(cfg["exclude_contigs"], list):
        cfg["exclude_contigs"] = []

    xchr = cfg.get("xchr")
    if xchr and isinstance(cfg["xchr"], str):
        cfg["xchr"] = xchr

    ychr = cfg.get("ychr")
    if ychr and isinstance(cfg["ychr"], str):
        cfg["ychr"] = ychr

    if cfg["include_contigs"]:
        log.debug(f"Trying to include these contigs in mosdepth: {', '.join(cfg['include_contigs'])}")
    if cfg["exclude_contigs"]:
        log.debug(f"Excluding these contigs from mosdepth: {', '.join(cfg['exclude_contigs'])}")
    if cfg.get("xchr"):
        log.debug(f'Using "{cfg["xchr"]}" as X chromosome name')
    if cfg.get("ychr"):
        log.debug(f'Using "{cfg["ychr"]}" as Y chromosome name')

    cutoff = cfg.get("perchrom_fraction_cutoff", 0.0)
    try:
        cutoff = float(cutoff)
    except ValueError:
        cutoff = 0.0
    if cutoff != 0.0:
        log.debug(f"Setting mosdepth coverage cutoff to display the contigs to {cutoff * 100.0}%")
    cfg["perchrom_fraction_cutoff"] = cutoff

    return cfg


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

    The MultiQC module plots coverage distributions from 2 kinds of outputs:

    - `{prefix}.mosdepth.region.dist.txt`
    - `{prefix}.mosdepth.global.dist.txt`

    Using "region" if exists, otherwise "global". Plotting 3 figures:

    - Proportion of bases in the reference genome with, at least, a given depth of coverage (cumulative coverage distribution).
    - Proportion of bases in the reference genome with a given depth of coverage (absolute coverage distribution).
    - Average coverage per contig/chromosome.

    Also plotting the percentage of the genome covered at a threshold in the General Stats section.
    The default thresholds are 1, 5, 10, 30, 50, which can be customised in the config as follows:

    ```yaml
    mosdepth_config:
      general_stats_coverage:
        - 10
        - 20
        - 40
        - 200
        - 30000
    ```

    You can also specify which columns would be hidden when the report loads (by default, all values are hidden except 30X):

    ```yaml
    general_stats_coverage_hidden:
      - 10
      - 20
      - 200
    ```

    For the per-contig coverage plot, you can include and exclude contigs based on name or pattern.

    For example, you could add the following to your MultiQC config file:

    ```yaml
    mosdepth_config:
      include_contigs:
        - "chr*"
      exclude_contigs:
        - "*_alt"
        - "*_decoy"
        - "*_random"
        - "chrUn*"
        - "HLA*"
        - "chrM"
        - "chrEBV"
    ```

    Note that exclusion superseeds inclusion for the contig filters.

    To additionally avoid cluttering the plot, mosdepth can exclude contigs with a low relative coverage.

    ```yaml
    mosdepth_config:
      # Should be a fraction, e.g. 0.001 (exclude contigs with 0.1% coverage of sum of
      # coverages across all contigs)
      perchrom_fraction_cutoff: 0.001
    ```

    If you want to see what is being excluded, you can set `show_excluded_debug_logs` to `True`:

    ```yaml
    mosdepth_config:
      show_excluded_debug_logs: True
    ```

    This will then print a debug log message (use `multiqc -v`) for each excluded contig.
    This is disabled by default as there can be very many in some cases.

    Besides the `{prefix}.mosdepth.global.dist.txt` and `{prefix}.mosdepth.region.dist.txt`
    files, the `{prefix}.mosdepth.summary.txt` file is used for the General Stats table.

    The module also plots an X/Y relative chromosome coverage per sample. By default, it finds chromosome named X/Y or chrX/chrY, but that can be customised:

    ```yaml
    mosdepth_config:
      # Name of the X and Y chromosomes. If not specified, MultiQC will search for
      # any chromosome names that look like x, y, chrx or chry (case-insensitive)
      xchr: myXchr
      ychr: myYchr
    ```
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="coverage",
            anchor="coverage",
            href="https://github.com/brentp/mosdepth",
            info="Fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing",
            doi="10.1093/bioinformatics/btx699",
        )

        self.cfg = read_config()
        genstats_by_sample: Dict[str, Dict[str, Union[int, float]]] = defaultdict(dict)  # mean coverage

        # Filter out any samples from --ignore-samples
        genstats_by_sample = defaultdict(dict, self.ignore_samples(genstats_by_sample))

        data_dicts_region = self.parse_cov_dist("region")
        data_dicts_region = [self.ignore_samples(d) for d in data_dicts_region]

        samples_region = set.union(*(set(d.keys()) for d in data_dicts_region))
        samples_found = samples_region
        if not samples_found:
            raise ModuleNoSamplesFound
        log.info(f"Found reports for {len(samples_found)} samples")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        descr_suf = ""
        if samples_region:
            descr_suf = ". Calculated across the target regions"

        if samples_region:
            data_dicts = data_dicts_region
            for d_region in data_dicts_region:
                d_region.update(d_region)
            (
                # cum_cov_dist_by_sample,
                # perchrom_avg_by_sample,
                #xy_cov_by_sample,
                extra_genstats_by_sample,
            ) = data_dicts

            # if cum_cov_dist_by_sample:
            #     xmax = 0
            #     for sample, cum_cov_by_x in cum_cov_dist_by_sample.items():
            #         for x, cumcov in cum_cov_by_x.items():
            #             if cumcov is not None and cumcov > 1:  # require >1% to prevent long flat tail
            #                 xmax = max(xmax, x)

            #     # Write data to file, sort columns numerically and convert to strings
            #     cumcov_dist_data_writeable = {
            #         sample: {str(k): v for k, v in sorted(cum_cov_by_x.items())}
            #         for sample, cum_cov_by_x in cum_cov_dist_by_sample.items()
            #     }
            #     self.write_data_file(cumcov_dist_data_writeable, "mosdepth_cumcov_dist")

            #     self.add_section(
            #         name="Cumulative coverage distribution",
            #         anchor="mosdepth-cumcoverage-dist",
            #         description=(
            #             f"Proportion of bases in the reference genome with, "
            #             f"at least, a given depth of coverage{descr_suf}"
            #         ),
            #         helptext=genome_fraction_helptext,
            #         plot=linegraph.plot(
            #             cum_cov_dist_by_sample,
            #             {
            #                 "id": "mosdepth-cumcoverage-dist-id",
            #                 "title": "Mosdepth: Cumulative coverage distribution",
            #                 "xlab": "Cumulative Coverage (X)",
            #                 "ylab": "% bases in genome/regions covered by at least X reads",
            #                 "ymin": 0,
            #                 "ymax": 100,
            #                 "xmin": 0,
            #                 "xmax": xmax,
            #                 "tt_label": "<b>{point.x}X</b>: {point.y:.2f}%",
            #                 "smooth_points": 500,
            #             },
            #         ),
            #     )

            #     # Write data to file, sort columns numerically and convert to strings
            #     cov_dist_data_writeable = {
            #         sample: {str(k): v for k, v in sorted(cum_cov_by_x.items())}
            #         for sample, cum_cov_by_x in cum_cov_dist_by_sample.items()
            #     }
            #     self.write_data_file(cov_dist_data_writeable, "mosdepth_cov_dist")

            # if perchrom_avg_by_sample:
            #     # Write data to file
            #     self.write_data_file(perchrom_avg_by_sample, "mosdepth_perchrom")

            #     num_contigs = max([len(x.keys()) for x in perchrom_avg_by_sample.values()])
            #     perchrom_plot: Union[Plot, str]
            #     if num_contigs > 1:
            #         perchrom_plot = linegraph.plot(
            #             perchrom_avg_by_sample,
            #             {
            #                 "id": "mosdepth-coverage-per-contig-multi",
            #                 "title": "Mosdepth: Coverage per contig",
            #                 "xlab": "Region",
            #                 "ylab": "Average Coverage",
            #                 "tt_decimals": 1,
            #                 "tt_suffix": "x",
            #                 "smooth_points": 500,
            #                 "logswitch": True,
            #                 "hide_empty": False,
            #                 "categories": True,
            #             },
            #         )
            #     else:
            #         perchrom_plot = bargraph.plot(
            #             perchrom_avg_by_sample,
            #             pconfig={
            #                 "id": "mosdepth-coverage-per-contig-single",
            #                 "title": "Mosdepth: Coverage per contig",
            #                 "xlab": "Sample",
            #                 "ylab": "Average Coverage",
            #                 "tt_suffix": "x",
            #                 "hide_empty": False,
            #                 "categories": True,
            #             },
            #         )

            #     self.add_section(
            #         name="Average coverage per contig",
            #         anchor="mosdepth-coverage-per-contig-section",
            #         description="Average coverage per contig or chromosome",
            #         plot=perchrom_plot,
            #     )

            # if xy_cov_by_sample:
            #     xy_keys = {
            #         "x": {"name": self.cfg.get("xchr") or "Chromosome X"},
            #         "y": {"name": self.cfg.get("xchr") or "Chromosome Y"},
            #     }
            #     pconfig = {
            #         "id": "mosdepth-xy-coverage-plot",
            #         "title": "Mosdepth: chrXY coverage",
            #         "ylab": "Coverage",
            #         "ysuffix": "X",
            #         "cpswitch_c_active": False,
            #     }
            #     self.add_section(
            #         name="XY coverage",
            #         anchor="mosdepth-xy-coverage",
            #         plot=bargraph.plot(
            #             {sname: {"x": x_cov, "y": y_cov} for sname, (x_cov, y_cov) in xy_cov_by_sample.items()},
            #             xy_keys,
            #             pconfig,
            #         ),
            #     )

            if extra_genstats_by_sample:
                update_dict(genstats_by_sample, extra_genstats_by_sample)

        # Adding hide buttons to report based of run names
        list_yaml_configs=["show_hide_buttons","show_hide_mode","show_hide_patterns,show_hide_color"]

        for yaml_header in list_yaml_configs:
            if not hasattr(config, yaml_header) or getattr(config,yaml_header) is None:
                setattr(config,yaml_header,{})
        list_show_hide_mode=[]
        list_show_hide_names=[]
        list_show_hide_color=[]
        more_less90={}

        for f in self.find_log_files(f"coverage/region_dist", filecontents=False, filehandles=True):
            s_name = self.clean_s_name(f["fn"], f)
            
            name=s_name.split("_",1)[0] if "_" in s_name else s_name
            if name not in list_show_hide_names:
                    list_show_hide_names.append(name)
                    list_show_hide_mode.append("show")
    
            if "_" not in s_name and s_name in genstats_by_sample:
                if name not in more_less90:
                    more_less90[name]=""
                if genstats_by_sample[s_name]['20_x_pc']<90:
                    more_less90[name]="#d9534f"

                if genstats_by_sample[s_name]['20_x_pc']>= 90 and more_less90[name]!="#d9534f":
                    more_less90[name]="#5cb85c"
        
        config.show_hide_buttons=list_show_hide_names
        config.show_hide_patterns=list_show_hide_names
        config.show_hide_mode=list_show_hide_mode
        for color in more_less90.values():
            list_show_hide_color.append(color)
        config.show_hide_color=list_show_hide_color

        # Generalstat_headers setup
        genstats_headers = {}
        threshs, hidden_threshs = config.get_cov_thresholds("mosdepth_config")
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
                "ceiling":90,
            }
        # Add mosdepth summary to General Stats
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
        self.general_stats_addcols(genstats_by_sample, genstats_headers)

    def parse_cov_dist(
        self, scope: str
    ) -> Tuple[
        Dict[str, Dict],
        Dict[str, Dict],
        Dict[str, Tuple[float, float]],
        Dict[str, Dict[str, Union[float, int, None]]],
    ]:
        """
        Two types of coverage distributions are parsed: global and region.

        {prefix}.mosdepth.global.dist.txt
        a distribution of proportion of bases covered at or above a given threshhold for each chromosome and genome-wide

        1       2       0.00
        1       1       0.00
        1       0       1.00
        total   2       0.00
        total   1       0.00
        total   0       1.00

        {prefix}.mosdepth.region.dist.txt (if --by is specified)
        same, but in regions

        1       2       0.01
        1       1       0.01
        1       0       1.00
        total   2       0.00
        total   1       0.00
        """

        # cumulative_pct_by_cov_by_sample: Dict[str, Dict] = defaultdict(dict)  # cumulative distribution
        bases_fraction_sum_per_contig_per_sample: Dict[str, Dict[str, float]] = defaultdict(
            dict
        )  # per chromosome average coverage
        #xy_cov_by_sample: Dict[str, Tuple[float, float]] = dict()
        genstats_by_sample: Dict[str, Dict[str, Union[float, int, None]]] = dict()

        threshs, hidden_threshs = config.get_cov_thresholds("mosdepth_config")

        # excluded_contigs = set()
        # included_contigs = set()
        # show_excluded_debug_logs = self.cfg.get("show_excluded_debug_logs") is True

        #Parse coverage distributions
        for f in self.find_log_files(f"coverage/{scope}_dist", filecontents=False, filehandles=True):
            s_name = self.clean_s_name(f["fn"], f)
                
            # if s_name in cumulative_pct_by_cov_by_sample:  # both region and global might exist, prioritizing region
            #     continue
            
            self.add_data_source(f, s_name=s_name, section="genome_results")

            bases_fraction_sum_per_contig: Dict[str, float] = defaultdict(float)
            cum_fraction_by_cov: Dict[int, float] = dict()

            for line in f["f"]:
                contig, cutoff_reads, bases_fraction = str(line).split("\t")
                if bases_fraction == "0.00\n":
                    continue

                # Parse cumulative coverage
                if contig == "total":
                    cum_fraction_by_cov[int(cutoff_reads)] = float(bases_fraction)

                # # Calculate per-contig coverage
                # else:
                #     # filter out contigs based on exclusion patterns
                #     if contig in excluded_contigs:
                #         continue

                #     if contig not in included_contigs:
                #         if any(fnmatch.fnmatch(contig, str(pattern)) for pattern in self.cfg["exclude_contigs"]):
                #             excluded_contigs.add(contig)
                #             if show_excluded_debug_logs:
                #                 log.debug(f"Skipping excluded contig '{contig}'")
                #             continue

                #         # filter out contigs based on inclusion patterns
                #         if len(self.cfg["include_contigs"]) > 0 and not any(
                #             fnmatch.fnmatch(contig, pattern) for pattern in self.cfg["include_contigs"]
                #         ):
                #             # Commented out since this could be many thousands of contigs!
                #             # log.debug(f"Skipping not included contig '{contig}'")
                #             continue

                #         included_contigs.add(contig)

                #     bases_fraction_sum_per_contig[contig] += float(bases_fraction)

            genstats_by_sample[s_name] = {}
            for k, v in genstats_cov_thresholds(cum_fraction_by_cov, threshs).items():
                genstats_by_sample[s_name][k] = v
                

            genstats_by_sample[s_name]["median_coverage"] = calc_median_coverage(cum_fraction_by_cov)

            # Downsampling the data to avoid carrying a lot for the line plot that would downsample anyway
            cum_fraction_by_cov = dict(smooth_array(list(cum_fraction_by_cov.items()), 500))
            # cumulative_pct_by_cov_by_sample[s_name] = {
            #     cutoff_reads: 100.0 * bases_fraction for cutoff_reads, bases_fraction in cum_fraction_by_cov.items()
            # }
            # bases_fraction_sum_per_contig_per_sample[s_name] = bases_fraction_sum_per_contig

        # Applying the contig coverage cutoff. First, count the total coverage for
        # every contig.
        total_cov_per_contig: Dict[str, float] = defaultdict(lambda: 0)
        total_cov = 0.0
        for s_name, bases_fraction_sum_per_contig in bases_fraction_sum_per_contig_per_sample.items():
            for contig, bases_fraction_sum in bases_fraction_sum_per_contig.items():
                total_cov_per_contig[contig] += bases_fraction_sum
                total_cov += bases_fraction_sum

        # Now, collecting the contigs that passed the cutoff.
        req_cov = float(total_cov) * self.cfg["perchrom_fraction_cutoff"]
        passing_contigs = set()
        for s_name, bases_fraction_sum_per_contig in bases_fraction_sum_per_contig_per_sample.items():
            for contig, bases_fraction_sum in bases_fraction_sum_per_contig.items():
                if float(total_cov_per_contig[contig]) > req_cov:
                    if contig not in passing_contigs:
                        passing_contigs.add(contig)

        rejected_contigs = set()
        filtered_perchrom_avg_data: Dict[str, Dict[str, float]] = defaultdict(dict)
        for s_name, bases_fraction_sum_per_contig in bases_fraction_sum_per_contig_per_sample.items():
            for contig, bases_fraction_sum in bases_fraction_sum_per_contig.items():
                if contig not in passing_contigs:
                    rejected_contigs.add(contig)
                else:
                    filtered_perchrom_avg_data[s_name][contig] = bases_fraction_sum
        bases_fraction_sum_per_contig_per_sample = filtered_perchrom_avg_data

        if rejected_contigs:
            if self.cfg.get("show_excluded_debug_logs") is True:
                log.debug(
                    f"Skipping {len(rejected_contigs)} contigs not passing the "
                    f"cutoff of {self.cfg['perchrom_fraction_cutoff']}% of "
                    f"{total_cov:.2f}x total coverage, which is {req_cov:.2f}x. "
                    f"Skipping contigs: {''.join(rejected_contigs)}"
                )

        # # Additionally, collect X and Y counts if we have them
        # for s_name, bases_fraction_sum_per_contig in bases_fraction_sum_per_contig_per_sample.items():
        #     x_cov: Optional[float] = None
        #     y_cov: Optional[float] = None
        #     for contig, bases_fraction_sum in bases_fraction_sum_per_contig.items():
        #         if self.cfg.get("xchr"):
        #             if str(self.cfg["xchr"]) == str(contig):
        #                 x_cov = bases_fraction_sum
        #         else:
        #             if contig.lower() == "x" or contig.lower() == "chrx":
        #                 x_cov = bases_fraction_sum
        #         if self.cfg.get("ychr"):
        #             if str(self.cfg["ychr"]) == str(contig):
        #                 y_cov = bases_fraction_sum
        #         else:
        #             if contig.lower() == "y" or contig.lower() == "chry":
        #                 y_cov = bases_fraction_sum
        #     # Only save these counts if we have both x and y
        #     if x_cov and y_cov:
        #         xy_cov_by_sample[s_name] = x_cov, y_cov

        # Correct per-contig average, since mosdepth reports cumulative coverage for at least
        # a certain value (see https://github.com/brentp/mosdepth#distribution-output).
        # For that reason, the 0 category (which is always 1) should not be included.
        for s_name in bases_fraction_sum_per_contig_per_sample:
            for contig in bases_fraction_sum_per_contig_per_sample[s_name]:
                bases_fraction_sum_per_contig_per_sample[s_name][contig] -= 1

        return (
            # cumulative_pct_by_cov_by_sample,
            # bases_fraction_sum_per_contig_per_sample,
            #xy_cov_by_sample,
            genstats_by_sample,
        )
