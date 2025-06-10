<p align="justify">
  <img src="docs/images/CMGG_logo.png" width="250" title="CMMG" alt="CMGG logo">
  
  <!-- Spacer: a transparent GIF or PNG -->
  <img src="https://via.placeholder.com/40x1/ffffff00?text=" width="5%" height="1" alt="" style="display:inline-block;">

  <picture>
    <source srcset="docs/images/MultiQC_logo_darkbg.png" media="(prefers-color-scheme: dark)">
    <source srcset="docs/images/MultiQC_logo.png" media="(prefers-color-scheme: light)">
    <img src="docs/images/MultiQC_logo.png" width="250" title="MultiQC" alt="MultiQC logo">
  </picture>
</p>

![Build status](https://github.com/ThomasSte01/MultiQC_CMGG/actions/workflows/build-docker.yaml/badge.svg) [![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/) 

**MultiQC_CMGG is a plugin for MultiQC containing customized modules and templates**

For more information about MultiQC, see [http://multiqc.info](http://multiqc.info)

## Description

The MultiQC_CMGG plugin add a custom template, used for internal QC reports and module customizations that can't be merged in core MultiQC.

## Installation

This plugin can be installed using the following methods

- using `pip`:

```bash
pip install --upgrade --force-reinstall git+https://github.com/CenterForMedicalGeneticsGhent/MultiQC_CMGG.git
```

- using `docker`:

```bash
docker build https://github.com/CenterForMedicalGeneticsGhent/MultiQC_CMGG.git -t multiqc_cmgg
```

## Usage

### Modules

| Name |Description| 
| ---- | --------- |
| sample_gender | Parse and merges ngs-bits [SampleGender](https://github.com/imgag/ngs-bits) output from xy, hetx, sry method into a sensible table. Based of [ngs-bits](https://docs.seqera.io/multiqc/modules/ngsbits) MultiQC module. To use module add `-m/--module sample_gender`. |
| coverage | Parse [Mosdepth](https://github.com/brentp/mosdepth) region files into tables based on genome panels. Based of [Mosdepth MultiQC module](https://docs.seqera.io/multiqc/modules/mosdepth) . To use module add `-m/--module coverage` and `-c/--config <path_to_multiqc_config_coverage.yaml>`|
| MSH2_hotspot_varcount |  Parses count files of a MSH2 hotspot variant from [nf-cmgg/report pipeline](https://github.com/nf-cmgg/report) to make a report in MultiQC. To use module add `-m/--module MSH2_hotspot_varcount` and `-c/--config <path_to_multiqc_config_msh2.yaml>` |
| msi_sensor_pro | Parses output from [msisensor-pro](https://github.com/xjtu-omics/msisensor-pro), generating an overview table from the summary file and a detailed table of MSI status for each locus from all results file. Also creates a bargraph visualizing the msisensor-pro scores. |

More thorough explaination on setting up modules, installation and integration into nextflow pipelines can be found [here](USAGE.md).

### Template

| Name | Description|
| ---- | ---------- |
| cmgg | CMGG specfic template with custom logo's and affiliate links. To enable this template, add the `-t/--template cmgg` option to the command line. |
