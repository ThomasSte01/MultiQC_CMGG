FROM python:latest

WORKDIR /home/multiqc_cmgg

COPY . .

# installing multiqc and multiqc_cmgg
RUN pip install multiqc
RUN pip install .


