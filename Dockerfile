FROM ubuntu:17.10
MAINTAINER Gareth Murphy <garethcmurphy@gmail.com>
RUN apt-get update && apt-get install -y --no-install-recommends  \
  ca-certificates=20161130+nmu1 \
  cmake \
  g++=4:6.3.0-2ubuntu2 \
  git=1:2.11.0-4 \
  libboost-all-dev=1.62.0.1  \
  libhdf5-serial-dev=1.10.0-patch1+docs-3 \
  make=4.1-9.1 \
  python3 \
  vim \
  && apt-get clean \
   && rm -rf /var/lib/apt/lists/*
COPY . /usr/src/myapp
WORKDIR /usr/src/myapp/build 
RUN cmake ..
RUN cmake --build .
RUN ctest -VV
