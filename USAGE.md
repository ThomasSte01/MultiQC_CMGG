# CenterForMedicalGeneticsGhent / multiqc_cmgg: Usage

## Modules

To ensure optimal performance and reduce report size, it's best practice to explicitly specify which modules should be used in your multiqc_config.yaml file. This prevents unnecessary modules from being loaded, which can slow down the process and bloat the final report.

This can be achieved by adding the following parameter to your multiqc_config.yaml file:

```yaml
run_modules:
  - coverage
  - MSH2_hotspot_varcount
  - sample_gender
  - ...
```

For each module a seperate multiqc_config.yaml template file is provided these can be found [here](docs/configs/).

### Coverage

The coverage module has certain parameters that should always be present, these being:

* no_violin: otherwise multiqc will convert this table into a violin plot when more then 500 rows of samples are present.
* coverage_config: to define which coverages should be calculated and hidden. 
* run_modules: because otherwise the build-in mosdepth module will be used
> **_Note:_** general_stats_coverage_hidden should not be empty, otherwise no coverage column will be shown

```yaml
custom_plot_config:
  general_stats_table:
    no_violin: true

coverage_config:
  general_stats_coverage: [20]
  general_stats_coverage_hidden: [30]

run_modules:
  - coverage
```

### MSH2_hotspot_varcount

The MSH2_hotspot_varcount module is designed to take results from the [nf-cmgg/report](https://github.com/nf-cmgg/report) pipeline which produces \<sample>.counts.txt files and parse these into a multiqc report.

This module requires the presence of the sanger_threshold parameter. This parameter defines the threshold at which samples are colored for failure or pass. 

```yaml
MSH2_hotspot_varcount_config:
  sanger_threshold : 20
```

### Sample_gender

The sample_gender modules does not require any parameters it will simply search for files containing *_sry.tsv , *_hetx.tsv ,and *_xy.tsv created by [ngs-bits SampleGender](https://github.com/imgag/ngs-bits/blob/master/doc/tools/SampleGender.md) and parse it into multiqc.

### Msi_sensor_pro



## Integrating into nextflow pipelines

To introduce this version of multiqc into a nextflow pipeline the nf-core module can be used.
First install the nf-core multiqc module with 
```bash 
nf-core modules install multiqc
```
Afterwards go to your modules.config file and add the following lines:

```groovy
process{
    withName: 'MULTIQC' {
        container ='quay.io/cmgg/multiqc_cmgg:0.0.2dev'
    }
}
```
This will change the docker image of multiqc to the available multiqc_cmgg image.

Next create or paste a multiqc_config.yaml file inside of the assets folder of your pipeline. 
Add the needed parameters inside of it as described in  the previous parts.
> **_Note:_** this file has to be called either multiqc_config.yaml or multiqc_config.yml

