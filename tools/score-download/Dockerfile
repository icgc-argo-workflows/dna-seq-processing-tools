FROM overture/score:1.6.2

RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository -y ppa:jonathonf/python-3.6 && \
    apt-get update && \
    apt-get install -y python3.6 && \
    ln -sf python3.6 /usr/bin/python3 && \
    ln -sf python3 /usr/bin/python


RUN chmod -R 777 /score-client

ENV PATH="/tools:/score-client/bin:${PATH}"

COPY *.py /tools/
