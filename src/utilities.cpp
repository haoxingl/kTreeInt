#include "utilities.h"
#include "config.h"

#include <cmath>
#include <chrono>

mpz_class generate_random_mpz(const mpz_class& min, const mpz_class& max) {
    static gmp_randclass rng(gmp_randinit_default);
    static bool seeded = false;
    if (!seeded) {
        unsigned long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        rng.seed(seed);
        seeded = true;
    }
    
    mpz_class range = max - min + 1;
    return min + rng.get_z_range(range);
}

mpf_class compute_p_denom(int m_exp, int k) {
    // Compute base exponent
    double base_exp = static_cast<double>(m_exp) / (log2(k) + 1);

    // Count divisions by 2
    int n_div = 0;
    while (base_exp > 32) {
        base_exp /= 2;
        n_div++;
    }

    // Compute 2^base_exp
    mpf_class result(0, config.mpf_precision);
    mpf_set_d(result.get_mpf_t(), exp2(base_exp));

    if (n_div > 0)
        mpf_pow_ui(result.get_mpf_t(), result.get_mpf_t(), 1 << n_div);

    return result;
}

mpz_class compute_m(int m_exp) {
    mpz_class m;
    mpz_ui_pow_ui(m.get_mpz_t(), 2, m_exp);
    return m;
}

mpf_class divide_mpz_by_mpf(const mpz_class& x, const mpf_class& y) {
    mpf_class result(x, config.mpf_precision);
    mpf_div(result.get_mpf_t(), result.get_mpf_t(), y.get_mpf_t());
    
    return result;
}

size_t compute_n(int m_exp, int k, double n_coefficient) {
    mpf_class n_f = compute_p_denom(m_exp, k);
    n_f *= n_coefficient;
    mpz_class n_z;
    mpz_set_f(n_z.get_mpz_t(), n_f.get_mpf_t());
    mpf_class diff = n_f - n_z;
    size_t n = n_z.get_ui();
    return diff < 0.5 ? n : n + 1;
}

bool is_valid_k(int m_exp, int k) {
    return k <= (1 << static_cast<size_t>((m_exp / log2(30)) - 1));
}

std::vector<double> generateSequence(int num_left, int num_right, double start_value, double interval) {
    std::vector<double> sequence(num_left + num_right + 1);
    sequence[0] = start_value;

    for (int i = 1; i <= num_left; ++i) {
        sequence[i] = sequence[i - 1] - interval;
    }
    sequence[num_left + 1] = start_value + interval;

    for (int i = num_left + 2; i <= num_left + num_right; ++i) {
        sequence[i] = sequence[i - 1] + interval;
    }

    return sequence;
}