<p float="left">
<img src="docs/images/CMGG_logo.png" width="250" title="CMMG">
<img src="docs/images/MultiQC_logo.png" width="250" title="MultiQC">
</p>

![Build status](https://github.com/CenterForMedicalGeneticsGhent/MultiQC_CMGG/actions/workflows/linux.yaml/badge.svg)

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

- using `setuptools`:

```bash
git clone https://github.com/CenterForMedicalGeneticsGhent/MultiQC_CMGG
cd MultiQC_CMGG
python setup.py install
```

## Usage

### Modules

| Name           | Description                                                                                                                                                                      |
| -------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| sample_gender | Parse and merges ngs-bits [SampleGender](https://github.com/imgag/ngs-bits) output from xy, hetx, sry method into a sensible table. Based of [ngs-bits](https://docs.seqera.io/multiqc/modules/ngsbits) MultiQC module. To use module add `-m/--module sample_gender`. |
| coverage | Parse Mosdepth [Mosdepth](https://github.com/brentp/mosdepth) region files into tables based on genome panels. Based of [Mosdepth](https://docs.seqera.io/multiqc/modules/mosdepth) MultiQC module. To use module add `-m/--module coverage`. |

### Templates

| Name | Description                                                                                                                                    |
| ---- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| cmgg | CMGG specfic template with custom logo's and affiliate links. To enable this template, add the `-t/--template cmgg` option to the command line. |
