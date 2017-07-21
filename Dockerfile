FROM ubuntu:17.10
MAINTAINER Gareth Murphy <garethcmurphy@gmail.com>
RUN apt-get update && apt-get install -y --no-install-recommends  \
  gnuplot \
  python3 \
  vim \
  ca-certificates \
  make \
  cmake    \
  g++ \
  git \
  libhdf5-serial-dev
COPY . /usr/src/myapp
WORKDIR /usr/src/myapp/build 
RUN cmake ..
RUN cmake --build .
RUN ctest -VV
