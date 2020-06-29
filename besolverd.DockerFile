FROM ubuntu:18.04

# Update the image to the latest packages
RUN apt update
RUN apt-get install -y libhts-dev;
RUN apt-get install -y python2.7;
RUN apt-get install -y tabix;
RUN apt-get install -y wget;
RUN apt-get install -y bedtools;
RUN apt-get install -y zip;

 
# INSTALL mosdepth
RUN mkdir /mosdepth
RUN wget 'https://github.com/brentp/mosdepth/releases/download/v0.2.9/mosdepth' -O mosdepth/mosdepth
RUN chmod a+x  mosdepth/mosdepth

# INSTALL rtg-tools
RUN wget 'https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.11/rtg-tools-3.11-linux-x64.zip'
RUN unzip rtg-tools-3.11-linux-x64.zip
RUN rm  rtg-tools-3.11-linux-x64.zip
COPY rtg.cfg /rtg-tools-3.11/

# Environment
ENV PATH="/rtg-tools-3.11:/mosdepth:${PATH}"

# SCRIPT DIR
RUN mkdir -p /besolverd/data
WORKDIR besolverd/data
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg37.bed
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg37chr.bed
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg37chr.vcf.gz
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg37chr.vcf.gz.tbi
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg37.vcf.gz
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg37.vcf.gz.tbi
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg38.bed
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg38.vcf.gz
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA12878_hg38.vcf.gz.tbi
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg37.bed
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg37chr.bed
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg37chr.vcf.gz
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg37chr.vcf.gz.tbi
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg37.vcf.gz
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg37.vcf.gz.tbi
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg38.bed
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg38.vcf.gz
RUN wget https://storage.googleapis.com/besolverd/giab_data/NA24385_hg38.vcf.gz.tbi

WORKDIR /besolverd
COPY besolverd.py .


