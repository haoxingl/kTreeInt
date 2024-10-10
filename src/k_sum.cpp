#include "k_sum.h"
#include "utilities.h"

// #include <random>
// #include <chrono>
// #include <stdexcept>
// #include <gmpxx.h>
#include <iostream>

bool DEBUG = 0;

std::vector<std::vector<mpz_class>> generate_k_sum_instance(const mpz_class& m, int k, size_t n) {
    if (k < 1 || n < 1) {
        throw std::invalid_argument("k and n must be positive integers");
    }

    std::vector<std::vector<mpz_class>> instance(k, std::vector<mpz_class>(n));

    if(DEBUG) {
        std::cout << "Generating k-SUM instance with m = " << m << ", k = " << k << ", n = " << n << std::endl;
    }

    for (int i = 0; i < k; ++i) {
        for (size_t j = 0; j < n; ++j) {
            instance[i][j] = generate_random_integer(m);
        }
    }

    return instance;
}

mpz_class generate_random_integer(const mpz_class& m) {
    return generate_random_mpz(-m / 2, m / 2);
}

bool check_sum_zero(const std::vector<mpz_class>& numbers) {
    mpz_class sum = 0;
    for (const auto& num : numbers) {
        sum += num;
    }
    return sum == 0;
}