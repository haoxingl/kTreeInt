#ifndef UTILITIES_H
#define UTILITIES_H

#include <gmpxx.h>
#include <vector>

// const Config config;

/**
 * @brief Generates a random mpz_class number within a given range.
 * 
 * @param min The minimum value (inclusive).
 * @param max The maximum value (inclusive).
 * @return mpz_class A random number between min and max.
 */
mpz_class generate_random_mpz(const mpz_class& min, const mpz_class& max);

/**
 * @brief Computes a power function based on given parameters.
 * 
 * This function computes 2^(m_exp / (log_2(k) + 1)) using a specific algorithm
 * to maintain precision with large exponents.
 *
 * @param m_exp The exponent numerator.
 * @param k The parameter used in the denominator of the exponent.
 * @return mpf_class The result as an mpf_class with 32-bit precision.
 */
mpf_class compute_p_denom(int m_exp, int k);

/**
 * @brief Computes m = 2^m_exp.
 * 
 * This function calculates m as 2 raised to the power of m_exp.
 *
 * @param m_exp The exponent.
 * @return mpz_class The result (2^m_exp) as an mpz_class.
 */
mpz_class compute_m(int m_exp);

/**
 * @brief Divides an mpz_class by an mpf_class.
 * 
 * This function divides x by y.
 *
 * @param x The mpz_class dividend.
 * @param y The mpf_class divisor.
 * @return mpf_class The result of x / y as an mpf_class.
 */
mpf_class divide_mpz_by_mpf(const mpz_class& x, const mpf_class& y);

/**
 * @brief Computes the value of n based on given parameters.
 * 
 * This function calculates n = n_coefficient * 2^(m_exp / (log_2(k) + 1)) * n_coefficient.
 *
 * @param m_exp The exponent of m.
 * @param k The number of lists.
 * @param n_coefficient The coefficient used to multiply the result.
 * @return size_t The result as an size_t.
 */
size_t compute_n(int m_exp, int k, double n_coefficient);

/**
 * @brief Checks if a given k value is valid for a given m_exp value.
 * 
 * This function checks if k is less than or equal to 2^(m_exp / log(30) - 1).
 *
 * @param m_exp The exponent of m.
 * @param k The number of lists.
 * @return bool True if k is valid, false otherwise.
 */
bool is_valid_k(int m_exp, int k);

std::vector<double> generateSequence(int num_left, int num_right, double start_value, double interval);

#endif // UTILITIES_H