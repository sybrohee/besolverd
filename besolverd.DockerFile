FROM ubuntu:bionic-20180426

# do not take into account the security or the updates of ubuntu
RUN mv /etc/apt/sources.list /etc/apt/sources.list.old
RUN cat /etc/apt/sources.list.old | grep -v '#' | grep -v 'updates' | grep -v 'security' > /etc/apt/sources.list

# Update the image to the latest packages
RUN apt update && apt install -y libhts-dev python tabix wget bedtools zip openjdk-8-jdk 


# INSTALL mosdepth
RUN wget 'https://github.com/brentp/mosdepth/releases/download/v0.2.9/mosdepth' && \
    chmod a+x  mosdepth && \
    mv mosdepth /usr/bin/

    
# PICARD installation
RUN mkdir -p /picard
RUN wget https://github.com/broadinstitute/picard/releases/download/2.23.4/picard.jar -O /picard/picard.jar

    
# INSTALL rtg-tools
RUN wget 'https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.11/rtg-tools-3.11-linux-x64.zip' && \
    unzip rtg-tools-3.11-linux-x64.zip && \
    rm rtg-tools-3.11-linux-x64.zip
#COPY rtg.cfg /rtg-tools-3.11/
RUN touch /rtg-tools-3.11/rtg.cfg

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




