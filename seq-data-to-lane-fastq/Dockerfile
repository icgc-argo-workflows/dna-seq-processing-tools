FROM ubuntu:20.04

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/dna-seq-processing-tools

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y && \
    apt-get install -y software-properties-common python3-pip python3-dev curl && \
    apt-get install -y libz-dev pkg-config libtool m4 autotools-dev automake libncurses5-dev libbz2-dev liblzma-dev

# install samtools
ARG SAMTOOLS_VERSION=1.15
RUN cd /tmp \
    && curl -sSL -o samtools-$SAMTOOLS_VERSION.tar.bz2 --retry 10 https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/samtools-$SAMTOOLS_VERSION.tar.bz2 \
    && bunzip2 -c samtools-$SAMTOOLS_VERSION.tar.bz2 |tar xf - \
    && cd samtools-$SAMTOOLS_VERSION \
    && ./configure --prefix=/usr/local \
    && make \
    && make install

ENV PATH="/tools:${PATH}"

COPY *.py /tools/

RUN groupadd -g 1000 ubuntu && \
    useradd -l -u 1000 -g ubuntu ubuntu && \
    install -d -m 0755 -o ubuntu -g ubuntu /home/ubuntu

ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
