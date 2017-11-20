#!/usr/bin/env bash
docker build --tag mhd .
docker run -v ${HOME}/test/mhdvanleer/build:/usr/src/myapp/build -it mhd
