#ifndef K_TREE_H
#define K_TREE_H

#include "observer.h"
#include "ResourceMonitor.h"

#include <gmpxx.h>
#include <vector>
#include <memory>


struct MergeElement {
    mpz_class value;
    std::vector<size_t> indices;
};

class KTree {
public:
    KTree(int k, int m_exp, mpf_class& p_denom);
    ~KTree() = default;

    std::vector<size_t> solve(std::vector<std::vector<mpz_class>>& lists);

    void add_observer(std::shared_ptr<Observer> observer);

private:

    int k;
    int m_exp;
    mpz_class m;
    mpf_class p_denom;
    std::vector<std::shared_ptr<Observer>> observers;
    
    size_t curr_level_size = 0;

    RuntimeMonitor runtimeMonitor;

    std::vector<MergeElement> merge(std::vector<MergeElement>& list1, std::vector<MergeElement>& list2, const mpf_class& tau);

    size_t binary_search(const std::vector<MergeElement>& sorted_list, size_t ind_start, size_t ind_end, const mpz_class& target, bool min_ind);

    void notify_observers(const std::string& event, const std::vector<size_t>& data);

    bool size_limit_exceeded(std::vector<MergeElement>& result, std::vector<MergeElement>& list1, std::vector<MergeElement>& list2);

    // extra statistics
    size_t count_zeros(const std::vector<MergeElement>& list);
};

#endif // K_TREE_H