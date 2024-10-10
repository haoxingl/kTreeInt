#ifndef OBSERVER_H
#define OBSERVER_H

#include <string>
#include <vector>


class Observer {
public:
    virtual ~Observer() = default;
    virtual void update(const std::string& event, const std::vector<size_t>& data) = 0;
};

class SuccessProbabilityObserver : public Observer {
public:
    SuccessProbabilityObserver();
    void update(const std::string& event, const std::vector<size_t>& data) override;
    double get_success_probability() const;
    int get_oversize_runs() const;
    bool size_within() const;

private:
    int total_runs;
    int successful_runs;
    int oversize_runs;
};

class ListSizeObserver : public Observer {
public:
    ListSizeObserver();
    void update(const std::string& event, const std::vector<size_t>& data) override;
    std::pair<double, double> get_total_size_avg_stddev() const;
    std::pair<double, double> get_max_level_size_avg_stddev() const;
    void print_all_total_size_and_dist_from_mean() const;

private:
    int total_runs;
    std::vector<size_t> total_sizes;
    std::vector<size_t> max_level_sizes;
    size_t current_level_size;
};


class ExtraLogObserver : public Observer {
public:
    ExtraLogObserver();
    void update(const std::string& event, const std::vector<size_t>& data) override;
    std::pair<double, double> get_zero_count_avg_stddev() const;
    double get_zero_count_second_moment() const;

private:
    int total_runs;
    std::vector<size_t> zero_counts;
};

#endif // OBSERVER_H