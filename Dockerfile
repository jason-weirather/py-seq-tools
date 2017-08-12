#py-seq-tools environment
FROM ubuntu:16.04
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y \
               python-pip \
    && apt-get autoremove \
    && apt-get clean

VOLUME /temp
VOLUME /root

RUN pip install --upgrade pip

RUN mkdir /source
COPY . /source/py-seq-tools

RUN cd /source/py-seq-tools && pip install -e .

ENV HOME /root
WORKDIR /root
