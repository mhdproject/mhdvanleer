# MHDVanleer

This is a serial version of a numerical code developed to solve the magnetohydrodynamics equations. It was developed primarily for jets in the interstellar medium. 

## Status

[![Build Status](https://travis-ci.org/garethcmurphy/mhdvanleer.svg?branch=master)](https://travis-ci.org/garethcmurphy/mhdvanleer)

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/9b7d2877bde0435d897f31e8c50497e6)](https://www.codacy.com/app/garethcmurphy/mhdvanleer?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=garethcmurphy/mhdvanleer&amp;utm_campaign=Badge_Grade)


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
hdf5, v1.8
```

### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be

```
docker build -t mhd .
```

And repeat

```
docker run -it mhd /bin/bash 
```

You can use 
```
python3 ../see.py 
```
to plot a file

## Running the tests

Google Test/ctest platform is used

```
ctest -VV
```



## Deployment

Add additional notes about how to deploy this on a live system

## Built With


## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Gareth Murphy** - *Initial work* - [garethcmurphy](https://github.com/garethcmurphy)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc
