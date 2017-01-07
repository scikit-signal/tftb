[![Documentation Status](https://readthedocs.org/projects/pytftb/badge/?version=master)](http://pytftb.readthedocs.org/en/master/?badge=master)
[![Build Status](https://travis-ci.org/scikit-signal/pytftb.svg)](https://travis-ci.org/scikit-signal/pytftb)
[![Coverage Status](https://coveralls.io/repos/scikit-signal/pytftb/badge.svg?branch=master&service=github)](https://coveralls.io/github/scikit-signal/pytftb?branch=master)
[![Code Health](https://landscape.io/github/scikit-signal/pytftb/master/landscape.svg?style=flat)](https://landscape.io/github/scikit-signal/pytftb/master)

pytftb
======

A Python implementation of the MATLAB Time-Frequency Toolbox by Auger, Flandrin, Goncalves and Lemoine (http://tftb.nongnu.org)


Documentation
-------------

Working draft of the documentation can be found [here](http://pytftb.rtfd.org).

Requirements
------------

1. Core requirements:
 * numpy
 * scipy
 * matplotlib
2. Optional requirements:
 * nose (for running tests)
 * sphinx (to build the documentation)


Installation
------------

1. Clone this repository with Git:

```bash
$ git clone https://github.com/scikit-signal/pytftb
```

2. Install the dependencies:

```bash
$ cd pytftb
$ pip install -r requirements.txt
```

3. Install PyTFTB:

```bash
$ python setup.py install
```
