#ifndef RESOURCE_MONITOR_H
#define RESOURCE_MONITOR_H

#include <chrono>
// #include <thread>
// #include <atomic>
// #include <functional>
#include "config.h"

class ResourceMonitor {
public:
    ResourceMonitor(double limit);
    virtual ~ResourceMonitor() = default;

    virtual double getCurrentUsage() const = 0;
    bool isLimitExceeded() const;

protected:
    virtual bool checkLimit() = 0;

    double limit;
};

class RuntimeMonitor : public ResourceMonitor {
public:
    RuntimeMonitor();
    RuntimeMonitor(double limit);
    double getCurrentUsage() const override;
    void reset();

protected:
    bool checkLimit() override;

private:
    std::chrono::steady_clock::time_point startTime;
};

#endif // RESOURCE_MONITOR_H