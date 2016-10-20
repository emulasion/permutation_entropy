#include <cmath>
#include <algorithm>
#include <map>

namespace permutation_entropy

{
    long int factorial(int num);

    long int recursive_permutation_rank(int length, int *vector, int *inverse);

    inline long int permutation_rank(int *vector, int length);

    double permutation_entropy_array_stats(double *, int, int);

    double permutation_entropy_dictionary_stats(double *, int, int);
}
