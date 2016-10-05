smartplus
=========



[![GitHub license](https://img.shields.io/badge/licence-GPL%203-blue.svg)](https://github.com/smartplus-team/smartplus/blob/master/LICENSE.txt)
[![Build Status](https://travis-ci.org/abetlen/smartplus.svg?branch=master)](https://travis-ci.org/abetlen/smartplus)

About
-----

smartplus is a free, open-source library for the simulation of heterogeneous materials. It is developed with the aim to be a high-quality scientific library to facilitate the analysis of the complex, non-linear composite material response and thus integrates several algorithms for the analysis of heterogeneous materials

smartplus is a C++ library with emphasis on speed and ease-of-use. Its principle focus is to provide tools to facilitate the implementation of up-to-date constitutive model for materials in Finite Element Analysis Packages. This is done by providing a C++ API to generate user material subroutine based on a library of functions. Also, SMART+ provides tools to analyse the behavior of material, considering loading at the material point level. Such tools include a thermomechanical solver, a software to predict effective properties of composites, and a built-in identification software (using a combined genetic-gradient based algorithm)

smartplus is mainly developed by contributors from the staff and students of Arts et Métiers ParisTech, that are members of the LEM3 laboratory. It is released under the GNU General Public License: GPL, version 3.

Documentation
--------------

Provider      | Status
--------      | ------
Read the Docs | [![Documentation Status](https://readthedocs.org/projects/smartplus/badge/?version=latest)](http://smartplus.readthedocs.io/en/latest)


Installation
------------

How to install SMART+ :

1 - Make sure you have Boost (1.60 at least) installed and Armadillo installed to use smartplus.

2 - Unzip the file in a source location, and rename it 'smartplus'.

3 - Go to such folder "smartplus"

4 - Execute the installation bash file : 

```bash
sh Install.sh
```

5 - Enjoy


How to use smartplus
--------------------

Several possibilities 

1 - Use the python wrapper for smartplus, simmit :

[![simmit](https://img.shields.io/badge/simmit-v0.9-blue.svg)](https://github.com/chemiskyy/simmit)

By doing, so, you can utilize smartplus in Ipython (jupyter) notebooks

2 - Use the executables provided by smartplus

For instance, a solver, an identification software, etc..

3 - Link smartplus with FEA Packages

An example is provided on how to use smartplus constitutive models as Umat libraries:

You can directly copy-paste "umat_single.o" (mechanical) or "umat_singleT.o" (thermomechanical) from 'pathtothefile'/smartplus/build/bin to your Abaqus work directory and use it like a classical Umat.

Example : 
```bash
abaqus job=mymodel.inp user=umat_single.o
```

4- Build your own projects using the SMART+ library

Link with -lsmartplus
Have fun :)

Authors
-------
* [Yves Chemisky](https://github.com/chemiskyy)
* [Kevin Bonnay](https://github.com/kbonnay)
