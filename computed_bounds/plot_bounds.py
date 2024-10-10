import argparse
import numpy as np
import matplotlib.pyplot as plt
from proof_functions import main_theorem, ProbBoundsAnaly, ShallueProbBounds, ShallueSizeBounds, JouxProbBounds, JouxSizeBounds
import warnings
# import matplotlib
# matplotlib.use('macosx')

PLOT_ANALY = True
type2_print_n = False
type2_plot_n = True

normal = "0123456789"
super_script = "⁰¹²³⁴⁵⁶⁷⁸⁹"
translator = str.maketrans(normal, super_script)

def binary_search_n(m_exp, k, prob_thres, p, zm, bound_type='ours', output_lb=False):
    left, right = 0.1, 1000.0
    
    if not output_lb:
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
            pass
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
                bound = a2 if not output_lb else a1
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

def plot_type2(m_exp, k_list, prob_thres, zm, bound_type, show, output_lb=False):
    plt.rcParams['axes.titlesize'] = 22  # Title font size
    plt.rcParams['axes.labelsize'] = 22  # X and Y label font size
    plt.rcParams['xtick.labelsize'] = 20  # X-axis tick label font size
    plt.rcParams['ytick.labelsize'] = 20  # Y-axis tick label font size
    plt.rcParams['legend.fontsize'] = 18  # Legend font size

    annotattion_fontsize = 16
    annotattion_ver_loc = 5
    figure_size = (8, 6)
    plot_markersize = 7.5
    
    y_lim_edge = 0.1

    size_lb_list = []
    size_ub_list = []
    valid_k_list = []
    n_list = []
    
    for k in k_list:
        log2_k = float(np.log2(k))
        p = 2.0 ** (-float(m_exp) / (log2_k + 1))
        
        c = binary_search_n(m_exp, k, prob_thres, p, zm, bound_type, output_lb=output_lb)
        if c is None:
            continue
        
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                ub, lb = calculate_bounds(k, c / p, m_exp, "size", zm, bound_type)
                size_ub_list.append(ub)
                size_lb_list.append(lb)
                valid_k_list.append(k)
                n_list.append(np.round(c/p))
        except (OverflowError, ValueError, ZeroDivisionError, RuntimeWarning):
            print(f"m = {m_exp}, skipping k={k} as it causes an exception when calculating size bounds")
            continue
        
    # check whether size_lb_list is empty
    if not size_lb_list:
        print("Error: No valid data points found")
        return
    
    bounds = size_ub_list if not output_lb else size_lb_list
            
    plt.figure(figsize=figure_size)
    x = np.log2(valid_k_list)
    
    plt.plot(x, np.log2(bounds), 'v' if not output_lb else '^', label='Complexity', color='blue', markersize=plot_markersize*1.5)
    
    for xi, yi in zip(x, np.log2(bounds)):
        log2_value = yi
        plt.annotate(f'{log2_value:.1f}', 
                     (xi, yi), 
                     xytext=(0, annotattion_ver_loc),
                     textcoords='offset points',
                     ha='center',
                     va='bottom',
                     fontsize=annotattion_fontsize)
    
    if type2_plot_n:
        plt.plot(x, np.log2(n_list), 'o', label='Input Size', color='green', markersize=plot_markersize)
        for xi, yi in zip(x, np.log2(n_list)):
            log2_value = yi
            plt.annotate(f'{log2_value:.1f}', 
                        (xi, yi), 
                        xytext=(0, -25),
                        textcoords='offset points',
                        ha='center',
                        va='bottom',
                        fontsize=annotattion_fontsize)
            
    # set ylim
    y_min = min(np.log2(n_list)) if type2_plot_n else min(np.log2(bounds))
    y_max = max(np.log2(bounds))
    plt.ylim(y_min - y_lim_edge * y_max, y_max + y_lim_edge * y_max)
    
    # Set x-axis to show only integer values
    min_log2k = int(np.floor(min(x)))
    max_log2k = int(np.ceil(max(x)))
    plt.xticks(range(min_log2k, max_log2k + 1))
    
    if type2_print_n:
        list_log2_k = range(min_log2k, max_log2k + 1)
        
        x_labels = [f'{x}\n[{np.log2(n):.2f}]' for x, n in zip(list_log2_k, n_list)]
        print(x_labels)
    
    plt.xlabel('log2(k)')
    plt.ylabel('Complexity (log2)')
    plt.title(f'{"Sufficient" if not output_lb else "Necessary"} Complexity: m = 2{str(int(m_exp)).translate(translator)}, Prob. = {prob_thres}')
    
    # transparent legend
    plt.legend(loc='upper right', framealpha=0.5)
    plt.grid(True)
    
    plt.tight_layout()
    zm_suffix = "_zm" if zm == 1 else ""
    bound_suffix = "_lb" if output_lb else "_ub"
    bound_type_suffix = f"_{bound_type}" if bound_type != 'ours' else ""
    plt.savefig(f'type2_m{int(m_exp)}_prob{prob_thres}{zm_suffix}{bound_suffix}{bound_type_suffix}.pdf', bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close()

def plot_type1(m_exp, k, zm, bound_type, show):
    
    plt.rcParams['axes.titlesize'] = 22  # Title font size
    plt.rcParams['axes.labelsize'] = 22  # X and Y label font size
    plt.rcParams['xtick.labelsize'] = 20  # X-axis tick label font size
    plt.rcParams['ytick.labelsize'] = 20  # Y-axis tick label font size
    plt.rcParams['legend.fontsize'] = 16  # Legend font size

    type2_flatten_threshold = 0.5e-3
    step_size = 0.001
    figure_size = (8, 6)
    plot_linewidth = 2.5
    NUM_XTICKS = 5
    
    try:
        log2_k = float(np.log2(k))
        p = 2.0 ** (-float(m_exp) / (log2_k + 1))
        
        # Find the starting point where upper bound becomes non-zero
        c_end = 0.1
        c_start = 0.05
        prev_lb = -1.0
        while True:
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings('error')
                    success_prob_ub, success_prob_lb = calculate_bounds(k, c_end / p, m_exp, "success_prob", zm, bound_type, remain_float=True)
                    if success_prob_ub > 0.01 and c_start < 0.1:
                        c_start = c_end
                    if c_end >= 1.0 and (success_prob_lb > 0.99 or success_prob_lb - prev_lb < type2_flatten_threshold):
                        break
                    c_end += step_size
                    prev_lb = success_prob_lb
            except (OverflowError, ValueError, ZeroDivisionError, RuntimeWarning):
                print(f"m={m_exp}, k={k}, Warning: exception occurred during finding valid range of c, c_start = {c_start}, c_end = {c_end}")
                if c_end > 0.1: c_end -= step_size
                break
                
        print(f"m={m_exp}, k={k}: c_start = {c_start}, c_end = {c_end}")

        c_values = np.linspace(c_start, c_end, 1000)
        n_values = c_values / p

        success_prob_ub_list = []
        success_prob_lb_list = []
        analy_prob_ub_list = []
        analy_prob_lb_list = []

        for n in n_values:
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings('error')
                    success_prob_ub, success_prob_lb = calculate_bounds(k, n, m_exp, "success_prob", zm, bound_type)
                    success_prob_ub_list.append(success_prob_ub)
                    success_prob_lb_list.append(success_prob_lb)
                    
                    if zm == 0 and PLOT_ANALY and bound_type == 'ours':
                        try:
                            with warnings.catch_warnings():
                                warnings.filterwarnings('error')
                                analy_prob_ub, analy_prob_lb = ProbBoundsAnaly(2**m_exp, k, n)
                                analy_prob_ub_list.append(analy_prob_ub)
                                analy_prob_lb_list.append(analy_prob_lb)
                        except (OverflowError, ValueError, ZeroDivisionError, RuntimeWarning):
                            pass
            except (OverflowError, ValueError, ZeroDivisionError, RuntimeWarning):
                print(f"Trying to compute the bounds: exception caught at c={c_values[len(success_prob_ub_list)]}")
                break

        if not success_prob_ub_list:
            print(f"Error: No valid data points for m={m_exp}, k={k}")
            return

        n_values = c_values / p

        # Plot and save the original graph
        plt.figure(figsize=figure_size)
        
        if zm == 0 and PLOT_ANALY and analy_prob_ub_list and bound_type == 'ours':
            # plt.figure(figsize=figure_size)
            plt.plot(n_values[:len(analy_prob_ub_list)], analy_prob_ub_list, 'b--', label='Thm 1 UB', linewidth=plot_linewidth)
            plt.plot(n_values[:len(analy_prob_lb_list)], analy_prob_lb_list, 'c--', label='Thm 1 LB', linewidth=plot_linewidth)
        plt.plot(n_values[:len(success_prob_ub_list)], success_prob_ub_list, 'g-', label='Thm 2 UB', linewidth=plot_linewidth)
        plt.plot(n_values[:len(success_prob_lb_list)], success_prob_lb_list, 'r-', label='Thm 2 LB', linewidth=plot_linewidth)
        
        x_min = n_values.min()
        x_max = n_values.max()
        x_ticks = np.linspace(x_min, x_max, NUM_XTICKS)
        c_ticks = x_ticks * p
        x_labels = [f'{np.log2(x):.3f}\n[{c:.3f}]' for x, c in zip(x_ticks, c_ticks)]

        plt.xticks(x_ticks, x_labels)

        plt.xlabel('Input List Size (log2) [c]')
        plt.ylabel('Succ. Prob.')
        plt.title(f'Succ. Prob. Bounds (m = 2{str(int(m_exp)).translate(translator)}, k = {k})')
        plt.legend(loc='upper left', framealpha=0.5)
        plt.grid(True)
        plt.tight_layout()
        zm_suffix = "_zm" if zm == 1 else ""
        bound_type_suffix = f"_{bound_type}" if bound_type != 'ours' else ""
        plt.savefig(f'type1_m{int(m_exp)}_k{k}{zm_suffix}{bound_type_suffix}.pdf', bbox_inches='tight')
        if show:
            plt.show()
        else:
            plt.close()

    except Exception as e:
        print(f"Error in plot_type1 for m={m_exp}, k={k}: {str(e)}")
        
def calculate_bounds(k, n, m_exp, mode, zm, bound_type, remain_float=False):
    m = 2**m_exp
    if bound_type == 'ours':
        return main_theorem(k, n, m_exp, mode, zm, remain_float=remain_float)
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

def main():
    parser = argparse.ArgumentParser(description="Generate plots using proof_functions.")
    parser.add_argument("-b", choices=['ours', 'shallue', 'joux'], default='ours', help="Bound type to use")
    parser.add_argument("-zm", type=int, choices=[0, 1], default=0, help="0 for original, 1 for Z_m version")
    parser.add_argument("-t", "--plot-type", choices=['type1', 'type2_ub', 'type2_lb'], required=True,
                        help="Plot type: 'type1', 'type2_ub' or 'type2_lb'")
    parser.add_argument("-m", type=int, required=True, help="m_exp value")
    parser.add_argument("-k", type=int, nargs='+', help="List of k values for type2 plot or single k value for type1 plot")
    parser.add_argument("-prob", type=float, help="Probability threshold for type2 plot")
    parser.add_argument("-show", type=int, choices=[0, 1], default=0, help="Show plot (1) or just save (0)")
    
    args = parser.parse_args()
    
    if args.plot_type == 'type2_ub':
        if not args.k or not args.prob:
            parser.error("type2 plot requires -k and -prob arguments")
        plot_type2(args.m, args.k, args.prob, args.zm, args.b, args.show)
    elif args.plot_type == 'type2_lb':
        if not args.k or not args.prob:
            parser.error("type2 plot requires -k and -prob arguments")
        plot_type2(args.m, args.k, args.prob, args.zm, args.b, args.show, True)
    elif args.plot_type == 'type1':
        if not args.k or len(args.k) != 1:
            parser.error("type1 plot requires a single -k argument")
        plot_type1(args.m, args.k[0], args.zm, args.b, args.show)
    else:
        parser.error("Invalid plot type")

if __name__ == "__main__":
    main()
