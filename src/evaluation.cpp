#include "evaluation.h"
#include "k_sum.h"
#include "k_tree.h"
#include "config.h"
#include "utilities.h"

#include <iostream>
#include <map>
#include <iomanip>
#include <ctime>
#include <filesystem>
#include <cmath>

double round_down_to_nearest_tenth(double value) {
    return std::floor(value * 100.0) / 100.0;
}

bool in_success_target(double target_prob, const std::vector<double>& success_target) {
    return std::find(success_target.begin(), success_target.end(), target_prob) != success_target.end();
}

Evaluation::Evaluation() 
    : runtimeMonitor() {
    // Generate result folder name based on current date and time
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S");
    result_folder = config.result_folder + "/" + ss.str();

    // Create result folder
    std::filesystem::create_directories(result_folder);

    // Open result file
    result_file.open(result_folder + "/results.csv");
    if (!result_file.is_open()) {
        throw std::runtime_error("Unable to open result file");
    }

    // Write metadata
    result_file << "# Evaluation started at: " << ss.str() << "\n";
    result_file << "# Global parameters:\n";
    result_file << "# TOTAL_SIZE_THRESHOLD: " << config.TOTAL_SIZE_THRESHOLD << "\n";
    result_file << "# MAX_RUNTIME_MINUTES: " << config.MAX_RUNTIME_MINUTES << "\n";
    result_file << "# MAX_RUNTIME_MINUTE_PERMERGE: " << config.MAX_RUNTIME_MINUTE_PERMERGE << "\n";
    result_file << "# total_runs_per_config: " << config.total_runs_per_config << "\n";
    // Add other relevant config parameters here
    result_file << "# mpf_precision: " << config.mpf_precision << "\n";

    // Write CSV header
    result_file << "m_exp,k,n_coefficient,success_probability,total_size_avg,total_size_stddev,max_level_size_avg,max_level_size_stddev,oversize_runs\n";

    // Open progress file
    progress_file.open(result_folder + "/progress.txt");
    if (!progress_file.is_open()) {
        throw std::runtime_error("Unable to open progress file");
    }

    if (config.EXTRA_LOG) {
        extra_logger.initialize_output();
        
        // Write metadata
        std::stringstream metadata;
        // metadata << "# Evaluation started at: " << extra_logger.generate_timestamp() << "\n";
        metadata << "# Global parameters:\n";
        metadata << "# TOTAL_SIZE_THRESHOLD: " << config.TOTAL_SIZE_THRESHOLD << "\n";
        metadata << "# MAX_RUNTIME_MINUTES: " << config.MAX_RUNTIME_MINUTES << "\n";
        metadata << "# MAX_RUNTIME_MINUTE_PERMERGE: " << config.MAX_RUNTIME_MINUTE_PERMERGE << "\n";
        metadata << "# total_runs_per_config: " << config.total_runs_per_config << "\n";
        metadata << "# mpf_precision: " << config.mpf_precision << "\n";
        extra_logger.write_metadata(metadata.str());
        
        extra_logger.write_header();
    }

    if (config.detailed_logging) {
        detailed_log_file.open(result_folder + "/" + config.detailed_log_file);
        if (!detailed_log_file.is_open()) {
            throw std::runtime_error("Unable to open detailed log file");
        }
    }
}

Evaluation::~Evaluation() {
    if (result_file.is_open()) {
        result_file.close();
    }
    if (progress_file.is_open()) {
        progress_file.close();
    }
    if (config.EXTRA_LOG) {
        extra_logger.close_files();
    }
    if (detailed_log_file.is_open()) {
        detailed_log_file.close();
    }
}

void Evaluation::write_result(const EvaluationResult& result) {
    result_file << result.m_exp << ","
                << result.k << ","
                << std::fixed << std::setprecision(16) << result.n_coefficient << ","
                << std::fixed << std::setprecision(3) << result.success_probability << ","
                << std::fixed << std::setprecision(3) << result.total_size_avg << ","
                << std::fixed << std::setprecision(3) << result.total_size_stddev << ","
                << std::fixed << std::setprecision(3) << result.max_level_size_avg << ","
                << std::fixed << std::setprecision(3) << result.max_level_size_stddev << ","
                << result.oversize_runs << "\n";
    result_file.flush();  // Periodically flush the file
    if (config.EXTRA_LOG) {
        // EvaluationResult latest_result = results.back();
        extra_logger.write_result(result.m_exp, result.k, result.n_coefficient, result.success_probability);
    }
}

void Evaluation::update_progress(const std::string& message) {
    progress_file << message << std::endl;
    progress_file.flush();
}

void Evaluation::run_evaluations() {
    std::cout << "Starting evaluation process..." << std::endl << std::endl;
    update_progress("Starting evaluation process...");

    for (int m_exp : config.m_exp_values) {
        for (int k : config.k_values) {
            if (!is_valid_k(m_exp, k)) continue;

            mpf_class p_denom = compute_p_denom(m_exp, k);
            bool permanent_skip_small_coef = false;

            for(double coef_step : config.n_coefficients_step) {
                std::vector<double> v_coef = generateSequence(9, 9, 1.0, coef_step);
                int success_n_coefficient_count = 0;
                bool skip_small_coef = permanent_skip_small_coef;
                for (double n_coefficient : v_coef) {
                    if (skip_small_coef && n_coefficient < 1.0) {
                        continue;
                    }

                    if (coef_step < 0.1 && n_coefficient == 1.0) continue;

                    bool oversize_flag = false;
                    bool overtime_flag = false;

                    try {
                        evaluate_single_configuration(m_exp, k, n_coefficient, p_denom);
                    } catch (const std::runtime_error& e) {
                        if (n_coefficient > 1.0 || (skip_small_coef && n_coefficient == 1.0)) {
                            std::cout << "Skipping configuration due to " << e.what() << std::endl;
                            update_progress("Skipping configuration due to " + std::string(e.what()));
                            break;  // Move to next k value
                        }
                        oversize_flag = e.what() == config.EVENT_OVERSIZE? true : oversize_flag;
                        overtime_flag = e.what() == config.EVENT_OVERTIME? true : overtime_flag;
                    }

                    if (!oversize_flag && !overtime_flag && !results.empty()) {
                        if (results.back().success_probability >= 0.001) {
                            write_result(results.back());
                            success_n_coefficient_count++;
                        }
                        if ((n_coefficient > 1.0 || (skip_small_coef && n_coefficient == 1.0)) && results.back().success_probability >= 0.999) {
                            std::cout << "Skipping evaluation early due to probability reached 1..." << std::endl;
                            update_progress("Skipping evaluation early due to probability reached 1...");
                            break;
                        }
                        if (n_coefficient <= 1.0 && results.back().success_probability < 0.001) {
                            if (n_coefficient == 1.0) permanent_skip_small_coef = true;
                            skip_small_coef = true; 
                        }
                    }
                }
                if (success_n_coefficient_count >= config.MIN_SUCCESS_C) break;
            }
        }
    }

    std::cout << "Evaluation process completed." << std::endl;
    update_progress("Evaluation process completed.");
}

void Evaluation::smart_run_evaluations() {
    log_detailed("Starting smart evaluation process...");

    for (int m_exp : config.m_exp_values) {
        for (int k : config.k_values) {

            std::unordered_map<double, EvaluationResult> run_history;

            mpf_class p_denom = compute_p_denom(m_exp, k);
            log_detailed("Evaluating m_exp: " + std::to_string(m_exp) + ", k: " + std::to_string(k));

            double next_target_threshold = 0.0;
            std::vector<double> success_target;

            for (double target_prob : config.success_prob_to_eval) {
                if (in_success_target(target_prob, success_target)) {
                    log_detailed("Skipping target probability " + std::to_string(target_prob) + " due to success target");
                    continue;
                }
                if (target_prob < next_target_threshold) {
                    log_detailed("Skipping target probability " + std::to_string(target_prob) + " due to threshold " + std::to_string(next_target_threshold));
                    continue;
                }

                log_detailed("\nTarget probability: " + std::to_string(target_prob));

                double n_left = config.bin_search_left;
                double n_right = config.bin_search_right;
                double n_mid = config.bin_search_start;

                // Check run history
                auto [key_left, key_right] = find_closest_keys(run_history, target_prob, target_prob == 1.0);
                log_detailed("Closest keys in run history: left = " + std::to_string(key_left) + ", right = " + std::to_string(key_right));

                // Adjust search range based on history
                if (key_left != -1.0) n_left = run_history[key_left].n_coefficient;
                if (key_right != 3.0) n_right = run_history[key_right].n_coefficient;
                log_detailed("Adjusted search range: left = " + std::to_string(n_left) + ", right = " + std::to_string(n_right));

                if (n_left != config.bin_search_left || n_right != config.bin_search_right) {
                    n_mid = (n_left + n_right) / 2;
                    log_detailed("Adjusted start point: n_mid = " + std::to_string(n_mid));
                }

                // Binary search
                int iteration = 0;
                bool oversize_flag = false;
                bool overtime_flag = false;
                while (true) {
                    iteration++;
                    log_detailed("Binary search iteration " + std::to_string(iteration) + ": n_mid = " + std::to_string(n_mid));

                    oversize_flag = false;
                    overtime_flag = false;
                    // EvaluationResult result;

                    try {
                        evaluate_single_configuration(m_exp, k, n_mid, p_denom);
                        log_detailed("Evaluation result: success_probability = " + std::to_string(results.back().success_probability));
                    } catch (const std::runtime_error& e) {
                        log_detailed("Caught runtime error: " + std::string(e.what()));
                        oversize_flag = (e.what() == config.EVENT_OVERSIZE);
                        overtime_flag = (e.what() == config.EVENT_OVERTIME);
                    }

                    if (!oversize_flag && !overtime_flag) {
                        double success_probability = results.back().success_probability;
                        if (success_probability >= 0.001) {
                            double rounded_prob = std::floor(success_probability * 1000.0) / 1000.0;
                            if (run_history.find(rounded_prob) == run_history.end()) {
                                run_history[rounded_prob] = results.back();
                                log_detailed("Added result to run history with key " + std::to_string(rounded_prob));
                            } else if (n_mid <= run_history[rounded_prob].n_coefficient) {
                                double prev_n_coefficient = run_history[rounded_prob].n_coefficient;
                                run_history[rounded_prob] = results.back();
                                log_detailed("Updated result in run history with key " + std::to_string(rounded_prob) + " with smaller n_coefficient " + std::to_string(n_mid) + " (prev: " + std::to_string(prev_n_coefficient) + ")");
                            }
                        }
                        
                        if ((target_prob != 1.0 && success_probability - target_prob <= config.bin_search_error_high && target_prob - success_probability <= config.bin_search_error_low) || 
                            (n_right - n_left) <= config.bin_search_min_interval) {
                            log_detailed("Search halted: last success_probability = " + std::to_string(success_probability));
                            if (config.DEBUG) {
                                log_detailed("n_left = " + std::to_string(n_left) + ", n_right = " + std::to_string(n_right));
                            }
                            for (double target_temp : config.success_prob_to_eval) {
                                if (target_temp == 1.0) continue;
                                if (std::find(success_target.begin(), success_target.end(), target_temp) == success_target.end()) {
                                    auto [key_left_temp, key_right_temp] = find_closest_keys(run_history, target_temp);
                                    if (key_left_temp != -1.0 || key_right_temp != 3.0) {
                                        if ((key_left_temp - target_temp <= config.bin_search_error_high && target_temp - key_left_temp <= config.bin_search_error_low)) {
                                            log_detailed("Adding success target " + std::to_string(target_temp));
                                            success_target.push_back(target_temp);
                                        } else if (key_right_temp - target_temp <= config.bin_search_error_high && target_temp - key_right_temp <= config.bin_search_error_low) {
                                            log_detailed("Adding success target " + std::to_string(target_temp));
                                            success_target.push_back(target_temp);
                                        }
                                    }
                                }
                            }
                            if ((n_right - n_left) <= config.bin_search_min_interval) {
                                auto [key_left, key_right] = find_closest_keys(run_history, target_prob);
                                log_detailed("Bin search halted due to interval limit for target " + std::to_string(target_prob) + ". Closest values in run history: left = " + std::to_string(key_left) + ", right = " + std::to_string(key_right));
                                next_target_threshold = round_down_to_nearest_tenth(key_right);
                                log_detailed("Next target threshold: " + std::to_string(next_target_threshold));
                            }
                            break;
                        }

                        if (success_probability < target_prob - config.bin_search_error_low) {
                            n_left = n_mid;
                            log_detailed("Adjusting left bound: new n_left = " + std::to_string(n_left));
                        } else {
                            n_right = n_mid;
                            log_detailed("Adjusting right bound: new n_right = " + std::to_string(n_right));
                        }
                    } else {
                        n_right = n_mid;
                        log_detailed("Error occurred. Adjusting right bound: new n_right = " + std::to_string(n_right));
                        if (run_history.find(2.0) == run_history.end() || 
                            n_mid < run_history[2.0].n_coefficient) {
                            EvaluationResult error_result = {m_exp, k, n_mid};
                            run_history[2.0] = error_result;
                            log_detailed("Updated error result in run history");
                        }
                        if (n_right - n_left <= config.bin_search_min_interval) {
                            log_detailed("Search halted due to interval limit");
                            if (config.DEBUG) {
                                log_detailed("n_left = " + std::to_string(n_left) + ", n_right = " + std::to_string(n_right));
                            }
                            break;
                        }
                    }

                    n_mid = (n_left + n_right) / 2;
                }
                if (oversize_flag || overtime_flag) {
                    log_detailed("Skipping the rest targets due to error flag");
                    if (config.DEBUG) {
                        log_detailed("Last target probability: " + std::to_string(target_prob));
                    }
                    break;
                } else {
                    double max_key = 0.0;
                    for (const auto& [key, value] : run_history) {
                        if (key <= 1.0 && key > max_key) max_key = key;
                    }
                    if (max_key < target_prob - config.bin_search_error_low) {
                        log_detailed("Skipping the rest targets due to low success probability");
                        if (config.DEBUG) {
                            log_detailed("Target probability: " + std::to_string(target_prob));
                            log_detailed("Success probability: " + std::to_string(results.back().success_probability));
                        }
                        break;
                    }
                }
            }
            for (double target_temp : config.success_prob_to_record) {
                auto [key_left_temp, key_right_temp] = find_closest_keys(run_history, target_temp);
                if (key_left_temp != -1.0 || key_right_temp != 3.0) {
                    if ((key_left_temp - target_temp <= config.bin_search_error_high && target_temp - key_left_temp <= config.bin_search_error_low) || (key_right_temp - target_temp <= config.bin_search_error_high && target_temp - key_right_temp <= config.bin_search_error_low)) {
                        double key_to_write = std::abs(key_left_temp - target_temp) <= std::abs(key_right_temp - target_temp) && (key_left_temp - target_temp <= config.bin_search_error_high && target_temp - key_left_temp <= config.bin_search_error_low) ? key_left_temp : key_right_temp;
                        EvaluationResult result_to_write = run_history[key_to_write];
                        log_detailed("Writing result for target probability " + std::to_string(target_temp) + " with key " + std::to_string(key_to_write));
                        write_result(result_to_write);
                    }
                }
            }
        }
    }

    log_detailed("Smart evaluation process completed.");
}

void Evaluation::evaluate_single_configuration(int m_exp, int k, double n_coefficient, mpf_class& p_denom) {
    std::cout << "Evaluating configuration: m_exp=" << m_exp 
    << ", k=" << k 
    << ", n_coefficient=" << std::fixed << std::setprecision(3) << n_coefficient 
    << std::endl;
    update_progress("Evaluating configuration: m_exp=" + std::to_string(m_exp) +
    ", k=" + std::to_string(k) + 
    ", n_coefficient=" + std::to_string(n_coefficient));

    mpz_class m = compute_m(m_exp);
    size_t n = compute_n(m_exp, k, n_coefficient);


    KTree ktree(k, m_exp, p_denom);
    std::shared_ptr<SuccessProbabilityObserver> success_observer = std::make_shared<SuccessProbabilityObserver>();
    std::shared_ptr<ListSizeObserver> list_size_observer = std::make_shared<ListSizeObserver>();
    ktree.add_observer(success_observer);
    ktree.add_observer(list_size_observer);

    // Extra log observer
    std::shared_ptr<ExtraLogObserver> extra_observer = std::make_shared<ExtraLogObserver>();
    ktree.add_observer(extra_observer);

    if (config.TIME_LIMIT) runtimeMonitor.reset();

    for (int run = 0; run < config.total_runs_per_config; ++run) {
        if (run % 100 == 0) {
            std::cout << "\rCompleted " << run << " out of " << config.total_runs_per_config << " trials" << std::flush;
        }

        auto instance = generate_k_sum_instance(m, k, n);
        ktree.solve(instance);
        std::vector<std::vector<mpz_class>>().swap(instance);

        if (config.TIME_LIMIT && runtimeMonitor.isLimitExceeded()) {
            std::cout << std::endl << "Stopping evaluation early due to runtime limit exceeded..." << std::endl;
            throw std::runtime_error(config.EVENT_OVERTIME);
        }
        if (!success_observer->size_within()) {
            std::cout << std::endl << "Stopping evaluation early due to too many oversize runs..." << std::endl;
            // break;
            throw std::runtime_error(config.EVENT_OVERSIZE);
        }
    }
    std::cout << std::endl;

    double success_probability = success_observer->get_success_probability();

    auto [total_size_avg, total_size_stddev] = list_size_observer->get_total_size_avg_stddev();
    auto [max_level_size_avg, max_level_size_stddev] = list_size_observer->get_max_level_size_avg_stddev();

    if (config.CHECK_STDDEV) {
        list_size_observer->print_all_total_size_and_dist_from_mean();
        std::cout << "Total size average: " << total_size_avg << ", Total size stddev: " << total_size_stddev << std::endl;
    }

    results.push_back({m_exp, k, n_coefficient, success_probability, total_size_avg, total_size_stddev, max_level_size_avg, max_level_size_stddev, success_observer->get_oversize_runs()});

    if (config.EXTRA_LOG) {
        auto [zero_count_avg, zero_count_stddev] = extra_observer->get_zero_count_avg_stddev();
        double second_moment = extra_observer->get_zero_count_second_moment();
        extra_logger.curr_zero_count = zero_count_avg;
        extra_logger.curr_second_moment = second_moment;
    }

    std::cout << "Finished trials for current configuration. Success probability: " << success_probability << std::endl << std::endl;
}

void Evaluation::log_detailed(const std::string& message) {
    if (config.DEBUG) std::cout << message << std::endl;
    if (config.detailed_logging) {
        detailed_log_file << message << std::endl;
        detailed_log_file.flush();
    }
}

std::pair<double, double> Evaluation::find_closest_keys(const std::unordered_map<double, EvaluationResult>& run_history, double target, bool is_one) {
    double key_left = -1.0;
    double key_right = 3.0;

    for (const auto& [key, value] : run_history) {
        if ((is_one? (key < target) : (key <= target)) && key > key_left) {
            key_left = key;
        }
        if (key >= target && key < key_right) {
            key_right = key;
        }
    }

    return {key_left, key_right};
}

