import argparse
import numpy as np
import matplotlib.pyplot as plt
from proof_functions import main_theorem, ShallueSizeBounds, JouxSizeBounds
import warnings

normal = "0123456789"
super_script = "⁰¹²³⁴⁵⁶⁷⁸⁹"
translator = str.maketrans(normal, super_script)

def binary_search_n(m_exp, k, prob_thres, p, zm, bound_type):
    left, right = 0.1, 1000.0
    
    # Check if the maximum n is still below the threshold or causes an exception
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            _, lb = calculate_bounds(k, right / p, m_exp, "success_prob", zm, bound_type)
        if lb < prob_thres:
            print(f"m = {m_exp}, skipping k={k} as very large c ({right}) doesn't reach the probability threshold (LB = {lb} < {prob_thres})")
            return None
    except (OverflowError, ValueError, ZeroDivisionError, RuntimeWarning):
        pass
    
    # Check if the minimum n causes an exception
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            _, _ = calculate_bounds(k, left / p, m_exp, "success_prob", zm, bound_type)
    except (OverflowError, ValueError, ZeroDivisionError, RuntimeWarning):
        print(f"m = {m_exp}, skipping k={k} as small c (p = {p}, c = {left}) already causes an overflow or invalid value")
        return None
    
    close_prob = -100.0
    max_prob = -100.0
    error_flag = False
    while right - left > 1e-3:
        error_flag = False
        mid = left + (right - left) / 2
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                a1, a2 = calculate_bounds(k, mid / p, m_exp, "success_prob", zm, bound_type)
                bound = a2
                if bound > max_prob:
                    max_prob = bound
                if bound >= prob_thres and abs(bound - prob_thres) < abs(close_prob - prob_thres):
                    close_prob = bound
            if bound < prob_thres:
                left = mid
            else:
                right = mid
        except (OverflowError, ValueError, ZeroDivisionError, RuntimeWarning):
            error_flag = True
            right = mid
    
    if not error_flag and close_prob >= prob_thres:
        return right
    
    print(f"m = {m_exp}, k={k}: Cannot find a c that required bound equals the probability threshold (closest value = {close_prob}, max value = {max_prob} vs {prob_thres}, c = {left}, {right})")
    return None

def calculate_bounds(k, n, m_exp, mode, zm, bound_type):
    m = 2**m_exp
    if bound_type == 'ours':
        return main_theorem(k, n, m_exp, mode, zm)
    elif bound_type == 'shallue':
        if mode == 'success_prob':
            return ShallueProbBounds(m, k, n)
        else:
            return ShallueSizeBounds(m, k, n)
    elif bound_type == 'joux':
        if mode == 'success_prob':
            return JouxProbBounds(m, k, n)
        else:
            return JouxSizeBounds(m, k, n)
    else:
        raise ValueError("Invalid bound type")

def plot_compare_type2(m_exp, k_list, prob_thres, zm, show):
    plt.rcParams['axes.titlesize'] = 18
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    plt.rcParams['legend.fontsize'] = 14

    annotation_fontsize = 18
    annotation_ver_loc = 5
    figure_size = (10, 8)
    plot_linewidth = 2.5
    plot_markersize = 7.5
    
    y_lim_edge = 0.2

    bounds = {
        'ours': [],
        'shallue': [],
        'joux': []
    }
    valid_k_list = []
    
    for k in k_list:
        log2_k = float(np.log2(k))
        p = 2.0 ** (-float(m_exp) / (log2_k + 1))
        
        c = binary_search_n(m_exp, k, prob_thres, p, zm, 'ours')
        if c is None:
            continue
        
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                for bound_type in bounds.keys():
                    ub, _ = calculate_bounds(k, c / p, m_exp, "size", zm, bound_type)
                    bounds[bound_type].append(ub)
                valid_k_list.append(k)
        except (OverflowError, ValueError, ZeroDivisionError, RuntimeWarning):
            print(f"m = {m_exp}, skipping k={k} as it causes an exception when calculating size bounds")
            continue
    
    if not any(bounds.values()):
        print("Error: No valid data points found")
        return
            
    plt.figure(figsize=figure_size)
    x = np.log2(valid_k_list)
    
    colors = {'ours': 'blue', 'shallue': 'green', 'joux': 'red'}
    markers = {'ours': 'v', 'shallue': 's', 'joux': 'o'}
    
    for bound_type, bound_list in bounds.items():
        label = ''
        curr_ver_loc = annotation_ver_loc
        if bound_type == 'ours':
            label = 'Theorem 2'
            curr_ver_loc = - 6*annotation_ver_loc
        if bound_type == 'shallue':
            label = '[Sha08]'
            curr_ver_loc = 5 * annotation_ver_loc
        if bound_type == 'joux':
            label = '[JKL24]'
            # curr_ver_loc = annotation_ver_loc
        plt.plot(x, np.log2(bound_list), markers[bound_type], 
                 #label=f'{bound_type.capitalize()} UB', 
                 label=label, 
                 color=colors[bound_type], 
                 markersize=plot_markersize,
                 linewidth=plot_linewidth)
        
        # Add annotations
        for xi, yi in zip(x, np.log2(bound_list)):
            log2_value = yi
            plt.annotate(f'{log2_value:.1f}',
                         (xi, yi),
                         xytext=(0, curr_ver_loc),
                         textcoords='offset points',
                         ha='center',
                         va='bottom',
                         fontsize=annotation_fontsize,
                         color=colors[bound_type])
    
    # Set y-axis limits
    y_min = min(min(np.log2(bound_list)) for bound_list in bounds.values() if bound_list)
    y_max = max(max(np.log2(bound_list)) for bound_list in bounds.values() if bound_list)
    plt.ylim(y_min - y_lim_edge * (y_max - y_min), y_max + y_lim_edge * (y_max - y_min))
    
    # Set x-axis to show only integer values
    min_log2k = int(np.floor(min(x)))
    max_log2k = int(np.ceil(max(x)))
    plt.xticks(range(min_log2k, max_log2k + 1))
    
    plt.xlabel('log2(k)')
    plt.ylabel('Complexity (log2)')
    plt.title(f'Sufficient Complexity Comparison: m = 2{str(int(m_exp)).translate(translator)}, Prob. = {prob_thres}')
    
    plt.legend(loc='upper right', framealpha=0.5)
    plt.grid(True)
    
    plt.tight_layout()
    zm_suffix = "_zm" if zm == 1 else ""
    plt.savefig(f'compare_type2_m{int(m_exp)}_prob{prob_thres}{zm_suffix}.pdf', bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close()

def main():
    parser = argparse.ArgumentParser(description="Generate comparison plots for type 2 bounds using different bound functions.")
    parser.add_argument("-zm", type=int, choices=[0, 1], default=0, help="0 for original, 1 for Z_m version")
    parser.add_argument("-m", type=int, required=True, help="m_exp value")
    parser.add_argument("-k", type=int, nargs='+', required=True, help="List of k values")
    parser.add_argument("-prob", type=float, required=True, help="Probability threshold")
    parser.add_argument("-show", type=int, choices=[0, 1], default=0, help="Show plot (1) or just save (0)")

    args = parser.parse_args()
    
    plot_compare_type2(args.m, args.k, args.prob, args.zm, args.show)

if __name__ == "__main__":
    main()
