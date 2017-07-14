FROM ubuntu:16.04
MAINTAINER Gareth Murphy <garethcmurphy@gmail.com>
RUN apt-get update && apt-get install -y --no-install-recommends  \
  gnuplot \
  python3 \
  vim \
  cmake \
  make \
  g++ \
  git \
  libhdf5-serial-dev \
  gcc
COPY . /usr/src/myapp
WORKDIR /usr/src/myapp/build 
RUN cmake ..
RUN cmake --build .
RUN ctest -VV
