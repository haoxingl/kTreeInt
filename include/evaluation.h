#ifndef EVALUATION_H
#define EVALUATION_H

#include "EvaluationLogger.h"
#include "ResourceMonitor.h"

#include <vector>
#include <gmpxx.h>
#include <fstream>
#include <unordered_map>


struct EvaluationResult {
    int m_exp;
    int k;
    double n_coefficient;
    double success_probability;
    double total_size_avg;
    double total_size_stddev;
    double max_level_size_avg;
    double max_level_size_stddev;
    int oversize_runs;
};

class Evaluation {
public:
    Evaluation();
    ~Evaluation();

    void run_evaluations();
    void smart_run_evaluations();

private:

    std::string result_folder;
    std::ofstream result_file;
    std::ofstream progress_file;
    std::vector<EvaluationResult> results;

    RuntimeMonitor runtimeMonitor;

    EvaluationLogger extra_logger;

    void write_result(const EvaluationResult& result);
    void update_progress(const std::string& message);
    void evaluate_single_configuration(int m_exp, int k, double n_coefficient, mpf_class& p_denom);

    std::ofstream detailed_log_file;
    void log_detailed(const std::string& message);

    std::pair<double, double> find_closest_keys(const std::unordered_map<double, EvaluationResult>& run_history, double target, bool is_one = false);
};

#endif // EVALUATION_H