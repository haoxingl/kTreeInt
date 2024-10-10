#include "EvaluationLogger.h"
#include "config.h"
#include <chrono>
#include <iomanip>
#include <sstream>
#include <filesystem>

EvaluationLogger::EvaluationLogger() {
    std::string timestamp = generate_timestamp();
    result_folder = config.result_folder + "/" + timestamp;
}

EvaluationLogger::EvaluationLogger(const std::string& base_folder) {
    std::string timestamp = generate_timestamp();
    result_folder = base_folder + "/" + timestamp;
}

EvaluationLogger::~EvaluationLogger() {
    close_files();
}

void EvaluationLogger::initialize_output() {
    std::filesystem::create_directories(result_folder);

    result_file.open(result_folder + "/results_extra.csv");
    if (!result_file.is_open()) {
        throw std::runtime_error("Unable to open result file");
    }
}

void EvaluationLogger::write_metadata(const std::string& metadata) {
    result_file << metadata;
    result_file.flush();
}

void EvaluationLogger::write_header() {
    result_file << "m_exp,k,n_coefficient,zero_count_avg,second_moment\n";
    result_file.flush();
}

void EvaluationLogger::write_result(const int m_exp, const int k, const double n_coefficient, double success_probability) {
    result_file << m_exp << ","
                << k << ","
                << std::fixed << std::setprecision(3) << n_coefficient << ","
                << std::fixed << std::setprecision(3) << success_probability << ","
                << std::fixed << std::setprecision(8) <<  curr_zero_count << ","
                << std::fixed << std::setprecision(8) <<  curr_second_moment
                << "\n";
    result_file.flush();
}

void EvaluationLogger::update_progress(const std::string& message) {
    progress_file << message << std::endl;
    progress_file.flush();
}

void EvaluationLogger::close_files() {
    if (result_file.is_open()) {
        result_file.close();
    }
    if (progress_file.is_open()) {
        progress_file.close();
    }
}

std::string EvaluationLogger::generate_timestamp() {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S");
    return ss.str();
}