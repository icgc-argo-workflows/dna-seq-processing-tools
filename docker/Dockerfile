FROM ubuntu:16.04

RUN apt-get update
RUN apt-get install -y python3-pip python3-dev
RUN apt-get install -y software-properties-common
RUN apt-get update
RUN apt-get install -y wget
RUN apt-get install -y samtools
RUN apt-get install -y bwa
RUN apt-get install -y git
RUN apt-get install -y libz-dev

RUN mkdir /tools

WORKDIR /tools
RUN wget https://github.com/biod/sambamba/releases/download/v0.6.9/sambamba-0.6.9-linux-static.gz
RUN gunzip sambamba-0.6.9-linux-static.gz
RUN mv sambamba-0.6.9-linux-static sambamba
RUN chmod +x /tools/sambamba
ENV PATH=/tools:$PATH
WORKDIR /

#Install Libmaus
#ARG LIBMAUS_VERSION=2.0.611-release-20190408112550
#WORKDIR /root
#RUN wget https://gitlab.com/german.tischler/libmaus2/-/archive/$LIBMAUS_VERSION/libmaus2-$LIBMAUS_VERSION.tar.gz
#RUN tar -zxvf libmaus2-$LIBMAUS_VERSION.tar.gz
#WORKDIR /root/libmaus2-$LIBMAUS_VERSION
#RUN ./configure
#RUN make
#WORKDIR /

#Install Biobambam
#ARG BIOBAMBAM_VERSION=2.0.95-release-20190320141403
#RUN apt-get install -y pkg-config
#WORKDIR /root
#RUN wget https://gitlab.com/german.tischler/biobambam2/-/archive/$BIOBAMBAM_VERSION/biobambam2-$BIOBAMBAM_VERSION.tar.gz
#RUN tar -zxvf biobambam2-$BIOBAMBAM_VERSION.tar.gz
#WORKDIR /root/biobambam2-$BIOBAMBAM_VERSION
#RUN ./configure --with-libmaus2=/root/libmaus2-$LIBMAUS_VERSION --prefix=${HOME}/biobambam2
#RUN make
#WORKDIR /

# Install OpenJDK-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get clean;

# Picard installation
RUN wget https://github.com/broadinstitute/picard/releases/download/2.19.0/picard.jar -O /tools/picard.jar
ENV PATH="/tools:${PATH}"

COPY . /app
WORKDIR /app