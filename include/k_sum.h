#ifndef K_SUM_H
#define K_SUM_H

#include <vector>
#include <gmpxx.h>

/**
 * @brief Generates a k-SUM input instance.
 * 
 * @param m The range of integers [-m/2, m/2] to draw from.
 * @param k The number of lists.
 * @param n The size of each list.
 * @return std::vector<std::vector<mpz_class>> A vector of k lists, each containing n integers.
 */
std::vector<std::vector<mpz_class>> generate_k_sum_instance(const mpz_class& m, int k, size_t n);

/**
 * @brief Generates a random integer in the range [-m/2, m/2].
 * 
 * @param m The upper bound of the range (exclusive).
 * @return mpz_class A random integer in the range [-m/2, m/2].
 */
mpz_class generate_random_integer(const mpz_class& m);

/**
 * @brief Checks if a given set of k integers sum to zero.
 * 
 * @param numbers A vector of k integers.
 * @return bool True if the sum is zero, false otherwise.
 */
bool check_sum_zero(const std::vector<mpz_class>& numbers);

#endif // K_SUM_H