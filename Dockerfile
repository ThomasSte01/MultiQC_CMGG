FROM quay.io/biocontainers/multiqc:1.29--pyhdfd78af_0

ADD . /src/multiqc_cmgg

RUN pip install --no-cache-dir -e /src/multiqc_cmgg

RUN multiqc --help
