FROM ubuntu:18.04

RUN apt-get update && \
    apt-get install -y python3-pip python3-dev software-properties-common wget && \
    apt-get install -y samtools bwa libz-dev pkg-config

RUN mkdir /tools

WORKDIR /tools
RUN wget https://github.com/biod/sambamba/releases/download/v0.6.9/sambamba-0.6.9-linux-static.gz
RUN gunzip sambamba-0.6.9-linux-static.gz
RUN mv sambamba-0.6.9-linux-static sambamba
RUN chmod +x /tools/sambamba
ENV PATH=/tools:$PATH

#Install libmaus2
LIBMAUS_VERSION=2.0.611-release-20190408112550
LIBMAUS_SOURCE=https://gitlab.com/german.tischler/libmaus2/-/archive/$LIBMAUS_VERSION/libmaus2-$LIBMAUS_VERSION.tar.gz
RUN cd /tmp \
    && wget LIBMAUS_SOURCE -O libmaus2-$LIBMAUS_VERSION.tar.gz \
    && tar -zxvf libmaus2-$LIBMAUS_VERSION.tar.gz \
    && cd libmaus2-$LIBMAUS_VERSION \
    && libtoolize \
    && aclocal \
    && autoreconf -i -f \
    && ./configure \
    && make \
    && make install

#Install biobambam2
BIOBAMBAM_VERSION=2.0.95-release-20190320141403
BIOBAMBAM_SOURCE=https://gitlab.com/german.tischler/biobambam2/-/archive/$BIOBAMBAM_VERSION/biobambam2-$BIOBAMBAM_VERSION.tar.gz
RUN cd /tmp \
    && wget BIOBAMBAM_SOURCE -O biobambam2-$BIOBAMBAM_VERSION.tar.gz \
    && tar -zxvf biobambam2-$BIOBAMBAM_VERSION.tar.gz \
    && cd biobambam2-$BIOBAMBAM_VERSION.tar.gz \
    && autoreconf -i -f \
    && ./configure --with-libmaus2=/usr/local/lib --prefix=/tools/biobambam2 \
    && make install

# Install OpenJDK-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get clean;

# Picard installation
RUN wget https://github.com/broadinstitute/picard/releases/download/2.19.0/picard.jar -O /tools/picard.jar

ENV LD_LIBRARY_PATH=/usr/local/lib
ENV PATH="/tools:/tools/biobambam2/bin:${PATH}"

COPY ./tools/* /tools/
CMD /bin/bash
