FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y gnuradio libboost-all-dev swig git

RUN git clone https://github.com/myriadrf/gr-limesdr

RUN apt-get install -y cmake

RUN apt-get install -y software-properties-common

RUN add-apt-repository -y ppa:myriadrf/drivers && apt-get update \
&& apt-get install -y limesuite liblimesuite-dev limesuite-udev limesuite-images \
soapysdr-tools soapysdr-module-lms7

RUN cd gr-limesdr && mkdir build && cd build && cmake .. && make && make install && ldconfig

RUN apt-get install -y python-pip && pip install numpy && pip install scipy

RUN apt-get install -y libcppunit-1.14-0 libfftw3-bin libvolk1-bin liblog4cpp5v5

RUN apt-get install -y python libliquid1d libliquid-dev

RUN git clone https://github.com/rpp0/gr-lora.git && cd gr-lora && mkdir build && cd build \
&& cmake .. && make && make install

RUN apt-get install -y xterm vim

WORKDIR /home

COPY ./flowgraphs /home
