from __future__ import unicode_literals
import itertools
import numpy as np
import sys
import math
import time

from pe_module_np import permutation_entropy as cpe

def log2(x):
    return np.log(x)/np.log(2)


def pype(time_series, m, delay=1):
    """Calculate the Permutation Entropy
    Args:
        time_series: Time series for analysis
        m: Order of permutation entropy
        delay: Time delay
    """
    n = len(time_series)
    permutations = np.array(list(itertools.permutations(range(m))))
    c = [0] * len(permutations)

    for i in xrange(n - delay * (m - 1)):
        # sorted_time_series =    np.sort(time_series[i:i+delay*m:delay], kind='quicksort')
        sorted_index_array = np.array(np.argsort(time_series[i:i + delay * m:delay], kind='quicksort'))
        for j in xrange(len(permutations)):
            if abs(permutations[j] - sorted_index_array).any() == 0:
                c[j] += 1

    c = [element for element in c if element != 0]
    p = np.array(c) / float(sum(c))
    pe = -sum(p * log2(p))
    return pe         


# generate some time series
# note that the time series must be a C contiguous array
time_series = np.ascontiguousarray(np.random.normal(0,1, 1000))

# embedding dimension
n = 6

print "C++ code, 'full_array_stat' method"
c_pe = cpe(time_series, n, "full_array_stat")
print "%.5f" % c_pe

print "C++ code, 'dict_array_stat' method"
c_pe = cpe(time_series, n, "dict_array_stat")
print "%.5f" % c_pe

print "pure python"
print "%.5f" % pype(time_series, n)        
        
        
        
        





