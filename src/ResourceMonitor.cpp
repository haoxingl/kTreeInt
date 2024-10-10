#include "ResourceMonitor.h"
#include <iostream>

ResourceMonitor::ResourceMonitor(double limit)
    : limit(limit) {}

bool ResourceMonitor::isLimitExceeded() const {
    return getCurrentUsage() > limit;
}

RuntimeMonitor::RuntimeMonitor()
    : ResourceMonitor(config.MAX_RUNTIME_MINUTES * 60) {
    startTime = std::chrono::steady_clock::now();
}

RuntimeMonitor::RuntimeMonitor(double limit)
    : ResourceMonitor(limit * 60) {
    startTime = std::chrono::steady_clock::now();
}

double RuntimeMonitor::getCurrentUsage() const {
    auto now = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(now - startTime).count();
}

void RuntimeMonitor::reset() {
    startTime = std::chrono::steady_clock::now();
}

bool RuntimeMonitor::checkLimit() {
    return getCurrentUsage() > limit;
}