#include "k_tree.h"
#include "utilities.h"
#include "config.h"

#include <cmath>
#include <stdexcept>
#include <iostream>

KTree::KTree(int k, int m_exp, mpf_class& p_denom) : k(k), m_exp(m_exp), p_denom(p_denom), runtimeMonitor(config.MAX_RUNTIME_MINUTE_PERMERGE) {
    m = compute_m(m_exp);
}

std::vector<size_t> KTree::solve(std::vector<std::vector<mpz_class>>& lists) {
    if (lists.size() != k) {
        throw std::invalid_argument("Number of lists must be equal to k");
    }

    mpf_class tau(m, config.mpf_precision);
    tau /= 2;

    notify_observers(config.EVENT_INIT_KLISTS, {k * lists[0].size()});

    std::vector<std::vector<MergeElement>> current_lists(k);
    for (int i = 0; i < k; ++i) {
        current_lists[i].reserve(lists[i].size());
        for (size_t j = 0; j < lists[i].size(); ++j) {
            current_lists[i].push_back({lists[i][j], {j}});
        }
        std::vector<mpz_class>().swap(lists[i]);
    }
    // swap lists to free memory
    std::vector<std::vector<mpz_class>>().swap(lists);

    for (int d = 1; d <= log2(k); ++d) {
        mpf_div(tau.get_mpf_t(), tau.get_mpf_t(), p_denom.get_mpf_t());
        
        std::vector<std::vector<MergeElement>> next_lists;
        curr_level_size = 0;

        for (size_t i = 0; i < current_lists.size(); i += 2) {
            std::vector<MergeElement> merged;
            try {
                merged = merge(current_lists[i], current_lists[i+1], tau);
            } catch (const std::exception& e) {
                if (e.what() == config.EVENT_MAX_SIZE_REACHED) {
                    std::vector<MergeElement>().swap(merged);
                    std::vector<std::vector<MergeElement>>().swap(next_lists);
                    std::vector<std::vector<MergeElement>>().swap(current_lists);
                    if (config.DEBUG) {
                        std::cerr << e.what() << '\n';
                    }
                    return {}; // Return empty vector to indicate failure
                }
            }
            

            std::vector<MergeElement>().swap(current_lists[i]);
            std::vector<MergeElement>().swap(current_lists[i+1]);
            
            if (merged.empty()) {
                notify_observers(config.EVENT_EMPTY_LIST, {static_cast<size_t>(d)});
                std::vector<MergeElement>().swap(merged);
                std::vector<std::vector<MergeElement>>().swap(next_lists);
                std::vector<std::vector<MergeElement>>().swap(current_lists);
                return {}; // Return empty vector to indicate failure
            }
            next_lists.push_back(merged);
            curr_level_size += merged.size();
            notify_observers(config.EVENT_MERGE, {merged.size()});
        }

        std::vector<std::vector<MergeElement>>().swap(current_lists);
        current_lists = next_lists;
        notify_observers(config.EVENT_LEVEL_COMPLETE,  {static_cast<size_t>(d)});
    }

    if (!current_lists[0].empty() && config.EXTRA_LOG) {
        size_t zeros = count_zeros(current_lists[0]);
        notify_observers(config.EVENT_EXTRA_LOG, {zeros});
    }

    // Check if 0 is in the final list
    auto it = std::find_if(current_lists[0].begin(), current_lists[0].end(), [](const MergeElement& elem) { return elem.value == 0; });
    if (it != current_lists[0].end()) {
        notify_observers(config.EVENT_SOLUTION_FOUND, {0});
        return it->indices;
    }

    notify_observers(config.EVENT_NO_SOLUTION, {});
    return {};
}

std::vector<MergeElement> KTree::merge(std::vector<MergeElement>& list1, std::vector<MergeElement>& list2, const mpf_class& tau) {

    if (config.SIZE_LIMIT) runtimeMonitor.reset();

    std::vector<MergeElement> result;

    if (config.DOUBLE_SORT) {
        std::sort(list1.begin(), list1.end(), [](const MergeElement& a, const MergeElement& b) { return a.value < b.value; });
    }
    std::sort(list2.begin(), list2.end(), [](const MergeElement& a, const MergeElement& b) { return a.value < b.value; });

    mpz_class floor_tau;
    mpz_set_f(floor_tau.get_mpz_t(), tau.get_mpf_t());

    size_t ptr2_left = list2.size() - 1;
    size_t ptr2_right = list2.size() - 1;

    for (const auto& elem1 : list1) {
        // Binary search for the last element c in list2 such that elem1.value + c.value <= floor_tau
        mpz_class target = floor_tau - elem1.value;
        size_t new_ptr2_right = binary_search(list2, 0, ptr2_right, target, false);

        // Binary search for the first element b in list2 such that elem1.value + b.value >= -floor_tau
        target = -floor_tau - elem1.value;
        size_t new_ptr2_left = binary_search(list2, 0, new_ptr2_right, target, true);

        if ((elem1.value + list2[new_ptr2_left].value) > floor_tau || elem1.value + list2[new_ptr2_right].value < -floor_tau) {
            continue;
        }

        if (config.SIZE_LIMIT && (curr_level_size + result.size() + new_ptr2_right - new_ptr2_left + 1) > config.TOTAL_SIZE_THRESHOLD) {
            size_limit_exceeded(result, list1, list2);
        }
        for (size_t i = new_ptr2_left; i <= new_ptr2_right; ++i) {
            MergeElement merged_elem;
            merged_elem.value = elem1.value + list2[i].value;
            merged_elem.indices = elem1.indices;
            merged_elem.indices.insert(merged_elem.indices.end(), list2[i].indices.begin(), list2[i].indices.end());
            result.push_back(merged_elem);

            if (config.SIZE_LIMIT && runtimeMonitor.isLimitExceeded()) {
                size_limit_exceeded(result, list1, list2);
            }
        }

        // Update pointers for the next iteration
        if (config.DOUBLE_SORT) {
            ptr2_left = new_ptr2_left;
            ptr2_right = new_ptr2_right;
        } else{
            size_t ptr2_left = list2.size() - 1;
            size_t ptr2_right = list2.size() - 1;
        }
    }

    return result;
}

size_t KTree::binary_search(const std::vector<MergeElement>& sorted_list, size_t ind_start, size_t ind_end, const mpz_class& target, bool min_ind) {
    // Check input validity
    if (ind_start < 0 || ind_end >= sorted_list.size() || ind_end < ind_start) {
        throw std::invalid_argument("Invalid input indices");
    }

    if (sorted_list[ind_start].value >= target) {
        return ind_start;
    }
    if (sorted_list[ind_end].value <= target) {
        return ind_end;
    }

    while (ind_end - ind_start > 1) {
        size_t ind_mid = ind_start + (ind_end - ind_start) / 2;
        
        if (sorted_list[ind_mid].value > target) {
            ind_end = ind_mid;
        } else if (sorted_list[ind_mid].value < target) {
            ind_start = ind_mid;
        } else {
            if (min_ind) {
                ind_end = ind_mid;
            } else {
                ind_start = ind_mid;
            }
        }
    }

    if (min_ind) {
        return (sorted_list[ind_start].value == target) ? ind_start : ind_end;
    } else {
        return (sorted_list[ind_end].value == target) ? ind_end : ind_start;
    }
}

void KTree::add_observer(std::shared_ptr<Observer> observer) {
    observers.push_back(observer);
}

void KTree::notify_observers(const std::string& event, const std::vector<size_t>& data) {
    for (auto& observer : observers) {
        observer->update(event, data);
    }
}

size_t KTree::count_zeros(const std::vector<MergeElement>& list) {
    size_t count = 0;
    for (const auto& elem : list) {
        if (elem.value == 0) {
            count++;
        }
    }
    return count;
}

bool KTree::size_limit_exceeded(std::vector<MergeElement>& result, std::vector<MergeElement>& list1, std::vector<MergeElement>& list2) {
    if (config.DEBUG) {
        std::cerr << "\nMerge: SIZE limit exceeded" << '\n';
        std::cerr << "Current level size: " << curr_level_size + result.size() << "/" << config.TOTAL_SIZE_THRESHOLD << '\n';
        std::cerr << "Runtime: " << runtimeMonitor.getCurrentUsage() << "/" << config.MAX_RUNTIME_MINUTE_PERMERGE * 60 << '\n';
    }
    notify_observers(config.EVENT_MAX_SIZE_REACHED, {result.size()});
    std::vector<MergeElement>().swap(result);
    std::vector<MergeElement>().swap(list1);
    std::vector<MergeElement>().swap(list2);
    throw std::runtime_error(config.EVENT_MAX_SIZE_REACHED);
}