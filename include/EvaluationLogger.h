#ifndef EVALUATION_LOGGER_H
#define EVALUATION_LOGGER_H

#include <string>
#include <fstream>
// #include "evaluation.h"  // For EvaluationResult struct

class EvaluationLogger {
public:
    EvaluationLogger();
    EvaluationLogger(const std::string& base_folder);
    ~EvaluationLogger();

    void initialize_output();
    void write_metadata(const std::string& metadata);
    void write_header();
    void write_result(const int m_exp, const int k, const double n_coefficient, double success_probability);
    void update_progress(const std::string& message);
    void close_files();

    double curr_zero_count = 0.0;
    double curr_second_moment = 0.0;

private:
    std::string result_folder;
    std::ofstream result_file;
    std::ofstream progress_file;

    std::string generate_timestamp();
};

#endif // EVALUATION_LOGGER_H