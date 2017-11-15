FROM ubuntu:17.10
MAINTAINER Gareth Murphy <garethcmurphy@gmail.com>
RUN apt-get update && apt-get install -y --no-install-recommends  \
  ca-certificates \
  cmake=3.9\*   \
  g++ \
  git \
  libboost-all-dev=1.62\*  \
  libhdf5-serial-dev=1.10\* \
  make=4.1\* \
  python3=3.6\*  \
  vim=2:8.0\* \
  && apt-get clean \
   && rm -rf /var/lib/apt/lists/*
COPY . /usr/src/myapp
WORKDIR /usr/src/myapp/build 
RUN cmake ..
RUN cmake --build .
RUN ctest -VV
