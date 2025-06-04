ARG MULTIQC_TAG=latest
FROM quay.io/biocontainers/multiqc:${MULTIQC_TAG}

ADD . /src/multiqc_cmgg

RUN pip install --no-cache-dir -e /src/multiqc_cmgg

RUN multiqc --help
