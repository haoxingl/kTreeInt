#include "observer.h"
#include "config.h"

#include <numeric>
#include <cmath>
// #include <chrono>
// #include <gmpxx.h>
#include <iostream>

// SuccessProbabilityObserver implementation
SuccessProbabilityObserver::SuccessProbabilityObserver() : total_runs(0), successful_runs(0), oversize_runs(0) {}

void SuccessProbabilityObserver::update(const std::string& event, const std::vector<size_t>& data) {
    if (event == config.EVENT_SOLUTION_FOUND) {
        successful_runs++;
    }
    if (event == config.EVENT_SOLUTION_FOUND || event == config.EVENT_NO_SOLUTION || event == config.EVENT_EMPTY_LIST) {
        total_runs++;
    }
    if (event == config.EVENT_MAX_SIZE_REACHED) {
        oversize_runs++;
    }
}

double SuccessProbabilityObserver::get_success_probability() const {
    return total_runs > 0 ? 
        static_cast<double>(successful_runs) / (total_runs) : 0.0;
}

int SuccessProbabilityObserver::get_oversize_runs() const {
    return oversize_runs;
}

bool SuccessProbabilityObserver::size_within() const {
    return oversize_runs < total_runs + 10;
}

ListSizeObserver::ListSizeObserver() : total_runs(0), total_sizes(std::vector<size_t>()), max_level_sizes(std::vector<size_t>()), current_level_size(0) {}

void ListSizeObserver::update(const std::string& event, const std::vector<size_t>& data) {
    if (event == config.EVENT_INIT_KLISTS){
        total_sizes.push_back(data[0]);
        max_level_sizes.push_back(data[0]);
        current_level_size = 0;
    } else if (event == config.EVENT_MERGE) {
        total_sizes.back() += data[0];
        current_level_size += data[0];
    } else if (event == config.EVENT_LEVEL_COMPLETE || event == config.EVENT_EMPTY_LIST) {
        if (current_level_size > max_level_sizes.back()) {
            max_level_sizes.back() = current_level_size;
        }
        current_level_size = 0;
        if (event == config.EVENT_EMPTY_LIST) {
            total_runs++;
        }
    } else if (event == config.EVENT_MAX_SIZE_REACHED) {
        if (total_sizes.size() > total_runs) {
            total_sizes.pop_back();
            max_level_sizes.pop_back();
        }
        current_level_size = 0;
    } else if (event == config.EVENT_SOLUTION_FOUND || event == config.EVENT_NO_SOLUTION) {
        total_runs++;
    }
}

std::pair<double, double> ListSizeObserver::get_total_size_avg_stddev() const {
    double avg = total_sizes.size() > 0 ? static_cast<double>(std::accumulate(total_sizes.begin(), total_sizes.end(), 0.0)) / total_sizes.size() : 0.0;
    double stddev = total_sizes.size() > 0 ? std::sqrt(std::accumulate(total_sizes.begin(), total_sizes.end(), 0.0, [avg](double acc, size_t x) { return acc + std::pow(x - avg, 2); }) / total_sizes.size()) : 0.0;
    return {avg, stddev};
}

std::pair<double, double> ListSizeObserver::get_max_level_size_avg_stddev() const {
    double avg = max_level_sizes.size() > 0 ? static_cast<double>(std::accumulate(max_level_sizes.begin(), max_level_sizes.end(), 0.0)) / max_level_sizes.size() : 0.0;
    double stddev = max_level_sizes.size() > 0 ? std::sqrt(std::accumulate(max_level_sizes.begin(), max_level_sizes.end(), 0.0, [avg](double acc, size_t x) { return acc + std::pow(x - avg, 2); }) / max_level_sizes.size()) : 0.0;
    return {avg, stddev};
}

void ListSizeObserver::print_all_total_size_and_dist_from_mean() const {
    double avg = get_total_size_avg_stddev().first;
    double stddev = get_total_size_avg_stddev().second;
    for (size_t i = 0; i < total_sizes.size(); i++) {
        std::cout << i+1 << ": Total size: " << total_sizes[i] << ", dist from mean: " << total_sizes[i] - avg << ", std dev: " << stddev << std::endl;
    }
}

ExtraLogObserver::ExtraLogObserver() : total_runs(0), zero_counts(std::vector<size_t>()) {}

void ExtraLogObserver::update(const std::string& event, const std::vector<size_t>& data) {
    if (event == config.EVENT_EXTRA_LOG) {
        zero_counts.push_back(data[0]);
    } else if (event == config.EVENT_SOLUTION_FOUND || event == config.EVENT_NO_SOLUTION || event == config.EVENT_EMPTY_LIST) {
        total_runs++;
    }
}

std::pair<double, double> ExtraLogObserver::get_zero_count_avg_stddev() const {
    double avg = zero_counts.size() > 0 ? static_cast<double>(std::accumulate(zero_counts.begin(), zero_counts.end(), 0.0)) / zero_counts.size() : 0.0;
    double stddev = zero_counts.size() > 0 ? std::sqrt(std::accumulate(zero_counts.begin(), zero_counts.end(), 0.0, [avg](double acc, size_t x) { return acc + std::pow(x - avg, 2); }) / zero_counts.size()) : 0.0;
    return {avg, stddev};
}

double ExtraLogObserver::get_zero_count_second_moment() const {
    return zero_counts.size() > 0 ? std::accumulate(zero_counts.begin(), zero_counts.end(), 0.0, [](double acc, size_t x) { return acc + x * x; }) / zero_counts.size() : 0.0;
}