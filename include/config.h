#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <string>

struct Config {
    bool DEBUG = 1;
    bool CHECK_STDDEV = 0;
    bool RAM_LIMIT = 0;
    bool TIME_LIMIT = 0;
    bool SIZE_LIMIT = 1;
    bool SMART_EVAL = 1;

    // K-Tree algorithm parameters
    size_t TOTAL_SIZE_THRESHOLD = 1ULL << 21; // GB /16
    size_t USED_RAM_THRESHOLD = 10 * 1024;
    size_t MERGE_RAM_CHECK_STEP = 1 * 1024;

    int total_runs_per_config = 1000;
    int mpf_precision = 17;
    std::string result_folder = "/home/hx/haoxing_k_tree/results";
    // Evaluation parameters
    
    double MAX_RUNTIME_MINUTES = 60.0;
    double MAX_RUNTIME_MINUTE_PERMERGE = 5.0/60;
    int MIN_SUCCESS_C = 10;
    // double CHECK_INTERVAL_SECONDS = 1.0;

    bool DOUBLE_SORT = false;
    bool EXTRA_LOG = false;
    
    std::string EVENT_SOLUTION_FOUND = "solution_found", 
                EVENT_NO_SOLUTION = "no_solution", 
                EVENT_EMPTY_LIST = "empty_list", 
                EVENT_MERGE = "merge", 
                EVENT_LEVEL_COMPLETE = "level_complete", 
                EVENT_INIT_KLISTS = "init_klists", 
                EVENT_MAX_SIZE_REACHED = "max_size_reached",
                EVENT_EXTRA_LOG = "extra_log",
                EVENT_OVERRAM = "overram",
                EVENT_OVERSIZE = "oversize",
                EVENT_OVERTIME = "overtime";

    // K-Sum problem parameters
    // std::vector<int> k_values = {1024, 512, 256, 128, 64, 32, 16, 8, 4};
    // std::vector<int> k_values = {1024, 512, 256, 128, 64, 32}; // m = 64 (1)
    // std::vector<int> k_values = {1024, 512, 256, 128, 64, 32}; // m = 96
    // std::vector<int> k_values = {1024}; // m = 128
    std::vector<int> k_values = {8, 16}; // m = 64 (2)

    // std::vector<int> m_exp_values = {64, 96, 128};
    std::vector<int> m_exp_values = {64};

    std::vector<double> n_coefficients_step = {0.1, 0.01, 0.001};

    // Smart Evaluation parameters
    // std::vector<double> success_prob_to_eval = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
    std::vector<double> success_prob_to_eval = {0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1.0};
    // std::vector<double> success_prob_to_eval = {0.5, 0.6, 0.7, 0.8, 0.9, 0.99};

    std::vector<double> success_prob_to_record = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,
    0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2,
    0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3,
    0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4,
    0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5,
    0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6,
    0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7,
    0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8,
    0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9,
    0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0};
    
    // double bin_search_left = 0.1;
    // double bin_search_left = 0.86; // m = 64 (1)
    // double bin_search_left = 0.85; // m = 96
    // double bin_search_left = 0.995; // m = 128
    double bin_search_left = 0.60; // m = 64 (2)
    // double bin_search_right = 1.9;
    // double bin_search_right = 1.12; // m = 64 (1)
    // double bin_search_right = 1.12; // m = 96
    // double bin_search_right = 1.005; // m = 128
    double bin_search_right = 1.46; // m = 64 (2)
    double bin_search_start = 1.0;
    double bin_search_error_high = 0.0099;
    double bin_search_error_low = 0.0;
    double bin_search_min_interval = 0.001;
    bool bin_search_always_add_valid_result = false;

    // Logging options
    bool detailed_logging = true;
    std::string detailed_log_file = "detailed_evaluation.log";

};

extern const Config config;

#endif // CONFIG_H