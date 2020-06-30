FROM debian:stretch-slim

# Update the image to the latest packages
RUN apt update && apt install -y libhts-dev python tabix wget bedtools zip


# INSTALL mosdepth
RUN wget 'https://github.com/brentp/mosdepth/releases/download/v0.2.9/mosdepth' && \
    chmod a+x  mosdepth && \
    mv mosdepth /usr/bin/

# INSTALL rtg-tools
RUN wget 'https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.11/rtg-tools-3.11-linux-x64.zip' && \
    unzip rtg-tools-3.11-linux-x64.zip && \
    rm rtg-tools-3.11-linux-x64.zip
#COPY rtg.cfg /rtg-tools-3.11/

# Environment
ENV PATH="/rtg-tools-3.11:${PATH}"

# SCRIPT DIR
RUN mkdir -p /besolverd/data && cd /besolverd/data  && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg37.bed && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg37chr.bed && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg37chr.vcf.gz && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg37chr.vcf.gz.tbi && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg37.vcf.gz && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg37.vcf.gz.tbi && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg38.bed && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg38.vcf.gz && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg38.vcf.gz.tbi && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg37.bed && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg37chr.bed && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg37chr.vcf.gz && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg37chr.vcf.gz.tbi && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg37.vcf.gz && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg37.vcf.gz.tbi && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg38.bed && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg38.vcf.gz && \
    wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg38.vcf.gz.tbi

WORKDIR /besolverd
ADD besolverd.py besolverd.py

ENTRYPOINT [ "python","besolverd.py" ]


