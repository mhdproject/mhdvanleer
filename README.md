# MHDVanleer

Simulating magnetised gas dynamics using shock capturing

## Status

[![Build Status](https://travis-ci.org/mhdproject/mhdvanleer.svg?branch=master)](https://travis-ci.org/mhdproject/mhdvanleer)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/9b7d2877bde0435d897f31e8c50497e6)](https://www.codacy.com/app/garethcmurphy/mhdvanleer?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=garethcmurphy/mhdvanleer&amp;utm_campaign=Badge_Grade)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them
Everything is installed with
```
Docker
```

Otherwise, dependencies are:

```
cmake (minimum 3.1)
g++, tested with v6.3, older versions may work
ctest
Google Test
hdf5, v1.10
```

### Installing

A step by step series of examples that tell you have to get a development env running

First build the Docker container and tag it:
```
docker build --tag mhd .
```
Then run the container

```
docker run -it mhd /bin/bash 
```

For visualisation, you can use 
```
python3 ../see.py 
```
to plot a file

## Running the tests

Google Test is automatically downloaded. The test runner is ctest.

```
ctest -VV
```



## Deployment

Add additional notes about how to deploy this on a live system

## Built With



[cmake](https://cmake.org/)

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/mhdproject/mhdvanleer/tags). 

## Authors

* **Gareth Murphy** - *Initial work* - [garethcmurphy](https://github.com/garethcmurphy)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc
