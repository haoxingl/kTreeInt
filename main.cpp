#include "evaluation.h"
// #include "config.h"
#include <iostream>
#include <chrono>

// const Config config;

int main() {

    std::cout << "Starting K-Tree Algorithm Evaluation" << std::endl;

    // Record start time
    auto start_time = std::chrono::high_resolution_clock::now();

    // Create and run evaluation
    Evaluation evaluator;
    if (config.SMART_EVAL) {
        evaluator.smart_run_evaluations();
    } else {
        evaluator.run_evaluations();
    }

    // Record end time
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

    std::cout << "Evaluation completed in " << duration.count() << " seconds." << std::endl;

    // You could add more summary statistics here if desired

    return 0;
}