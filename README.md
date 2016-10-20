# permutation_entropy

Small C++ library for fast evaluation of the permutation entropy of a time series (see [this paper](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.88.174102) for additional details). A [NumPy](http://www.numpy.org/) wrapper is included.

For a quick test:

compile the python module
```
python setup.py build_ext --inplace
```
then lauch the python test script (the NumPy library is needed):
```
python test.py
```
