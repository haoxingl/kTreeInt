import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os
import argparse

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

# import all functions from proof_functions.py
from proof_functions import main_theorem, chernoff_bounds

normal = "0123456789"
super_script = "⁰¹²³⁴⁵⁶⁷⁸⁹"
translator = str.maketrans(normal, super_script)

def plot_type1(group, m_exp, k, pdf=None, output_path=None):
    plt.rcParams['axes.titlesize'] = 22  # Title font size
    plt.rcParams['axes.labelsize'] = 22  # X and Y label font size
    plt.rcParams['xtick.labelsize'] = 20  # X-axis tick label font size
    plt.rcParams['ytick.labelsize'] = 20  # Y-axis tick label font size
    plt.rcParams['legend.fontsize'] = 16  # Legend font size

    figure_size = (8, 6)
    plot_linewidth = 2.5
    plot_markersize = 7.5
    NUM_XTICKS = 5
    
    if len(group) < 5:
        print(f"Skipping group with m_exp={m_exp}, k={k} due to insufficient data points")
        return  # Skip this group
    
    # Ensure m_exp and k are native Python types
    m_exp = float(m_exp)
    k = int(k)
    log2_k = float(np.log2(k))

    # Compute p using the new definition
    p = 2.0 ** ( - float(m_exp) / (float(log2_k) + 1) )

    # Compute c_values and n_values
    c_values = np.linspace(group['n_coefficient'].min() * 0.999, group['n_coefficient'].max() * 1.001, 1000)

    # Sort and filter the group data
    filtered_group = group.sort_values('n_coefficient')
    max_prob = float('-inf')
    filtered_indices = []
    for idx, row in filtered_group.iterrows():
        if row['success_probability'] > max_prob:
            filtered_indices.append(idx)
            max_prob = row['success_probability']
    filtered_group = filtered_group.loc[filtered_indices]

    if len(filtered_group) < 5:
        print(f"Skipping group with m_exp={m_exp}, k={k} due to insufficient filtered data points")
        return  # Skip this group

    # Calculate Chernoff bounds
    n_trials = 1000  # number of trials for each point
    confidence = 0.99  # 99% confidence interval
    chernoff_lower = []
    chernoff_upper = []
    for prob in filtered_group['success_probability']:
        lower, upper = chernoff_bounds(prob, n_trials, confidence)
        chernoff_lower.append(lower)
        chernoff_upper.append(upper)

    # Calculate theoretical upper and lower bounds using the new functions
    success_prob_ub_list = []
    success_prob_lb_list = []
    for n_coeff in c_values:
        n = n_coeff / p
        success_prob_ub, success_prob_lb = main_theorem(k, n, m_exp, "success_prob")
        success_prob_ub_list.append(success_prob_ub)
        success_prob_lb_list.append(success_prob_lb)

    plt.figure(figsize=figure_size)
    plt.plot(filtered_group['n_coefficient']/p, filtered_group['success_probability'], 'bo', label='Empirical', markersize=plot_markersize)

    plt.plot(c_values/p, success_prob_ub_list, 'g--', label='UB', linewidth=plot_linewidth)
    plt.plot(c_values/p, success_prob_lb_list, 'r-.', label='LB', linewidth=plot_linewidth)

    # Add Chernoff bounds
    lower_error = filtered_group['success_probability'] - chernoff_lower
    upper_error = np.array(chernoff_upper) - filtered_group['success_probability']

    # Plot the error bars
    plt.errorbar(filtered_group['n_coefficient']/p, filtered_group['success_probability'], 
                yerr=[lower_error, upper_error], fmt='none', ecolor='blue', 
                capsize=5, label='99% Conf.', alpha=0.5)

    # Limit x-axis based on n_coefficient
    x_min = group['n_coefficient'].min() * 0.999 / p
    x_max = group['n_coefficient'].max() * 1.001 / p
    
    # Create custom x-ticks with c_values in brackets
    num_xticks = NUM_XTICKS  # You can adjust this number
    x_ticks = np.linspace(x_min, x_max, num_xticks)
    c_ticks = x_ticks * p
    x_labels = [f'{np.log2(x):.3f}\n[{c:.3f}]' for x, c in zip(x_ticks, c_ticks)]

    plt.xticks(x_ticks, x_labels)

    plt.xlabel('Input List Size (log2) [c]')
    plt.ylabel('Succ. Prob.')
    plt.title(f'Succ. Prob. vs Input Size (m = 2{str(int(m_exp)).translate(translator)}, k = {k})')
    plt.legend(loc='upper left')
    plt.grid(True)
    
    plt.subplots_adjust(bottom=0.2)
    
    if pdf:
        pdf.savefig(bbox_inches='tight')
    elif output_path:
        plt.savefig(output_path, bbox_inches='tight')
    plt.close()

def plot_type2(group, m_exp, success_probability, pdf=None, output_path=None):
    plt.rcParams['axes.titlesize'] = 22  # Title font size
    plt.rcParams['axes.labelsize'] = 22  # X and Y label font size
    plt.rcParams['xtick.labelsize'] = 20  # X-axis tick label font size
    plt.rcParams['ytick.labelsize'] = 20  # Y-axis tick label font size
    plt.rcParams['legend.fontsize'] = 16  # Legend font size

    annotattion_fontsize = 16
    annotattion_ver_loc = 3
    figure_size = (8, 6)
    plot_linewidth = 2.5
    plot_markersize = 7.5
    NUM_YTICKS = 10
    
    k_start = 1
    
    def get_first_above_threshold(group):
        above_threshold = group[group['success_probability'] >= success_probability]
        if not above_threshold.empty:
            return above_threshold.iloc[0]
        return None

    filtered_group = group.groupby('k').apply(get_first_above_threshold)
    filtered_group = filtered_group.dropna()

    if filtered_group.empty:
        print(f"Skipping plot_type2 for m_exp={m_exp} and success_probability={success_probability} due to insufficient data")
        return

    fig, ax = plt.subplots(figsize=figure_size)

    x = np.log2(filtered_group.index)
    
    m_exp = float(m_exp)

    # Calculate p for each k
    p_values = filtered_group.apply(lambda row: 2.0 ** ( - float(m_exp) / (float(np.log2(row.name)) + 1) ), axis=1)
    input_size = filtered_group['n_coefficient'] / p_values
    
    ax.plot(x, filtered_group['total_size_avg'], 'o', label='Time', color='green',markersize=plot_markersize)
    for xi, yi in zip(x, filtered_group['total_size_avg']):
        log2_value = np.log2(yi)
        ax.annotate(f'{log2_value:.1f}', 
                    (xi, yi), 
                    xytext=(0, annotattion_ver_loc),
                    textcoords='offset points',
                    ha='center',
                    va='bottom',
                    fontsize=annotattion_fontsize)

    # Calculate and plot size bounds using the new functions
    size_ub_values = []
    size_lb_values = []
    for idx, row in filtered_group.iterrows():
        k = row.name
        n = row['n_coefficient'] / p_values[idx]
        size_ub, size_lb = main_theorem(k, n, m_exp, "size")
        size_ub_values.append(size_ub)
        size_lb_values.append(size_lb)

    ax.plot(x, size_ub_values, '--', label='Time UB', color='darkgreen', linewidth=plot_linewidth)
    ax.plot(x, size_lb_values, '-.', label='Time LB', color='lightgreen', linewidth=plot_linewidth)
    
    ax.plot(x, filtered_group['max_level_size_avg'], 'o', label='Space', color='red',markersize=plot_markersize)
    for xi, yi in zip(x, filtered_group['max_level_size_avg']):
        log2_value = np.log2(yi)
        ax.annotate(f'{log2_value:.1f}', 
                    (xi, yi), 
                    xytext=(0, annotattion_ver_loc),
                    textcoords='offset points',
                    ha='center',
                    va='bottom',
                    fontsize=annotattion_fontsize)

    # Plot lines instead of bars
    ax.plot(x, input_size, 'o', label='Input Size', color='blue', markersize=plot_markersize)
    for xi, yi in zip(x, input_size):
        log2_value = np.log2(yi)
        ax.annotate(f'{log2_value:.1f}', 
                    (xi, yi), 
                    xytext=(0, annotattion_ver_loc),
                    textcoords='offset points',
                    ha='center',
                    va='bottom',
                    fontsize=annotattion_fontsize)

    # Add error bars for total_size_avg and max_level_size_avg
    ax.errorbar(x, filtered_group['total_size_avg'], yerr=filtered_group['total_size_stddev'], fmt='o',capsize=5, ecolor='green', alpha=0.5)
    ax.errorbar(x, filtered_group['max_level_size_avg'], yerr=filtered_group['max_level_size_stddev'], fmt='o', capsize=5, ecolor='red', alpha=0.5)

    # Set y-axis limit
    ax.set_ylim(-0.05 * filtered_group['total_size_avg'].max(), 1.1 * filtered_group['total_size_avg'].max())

    # Set x-axis to show only integer values
    min_log2k = int(np.floor(min(x)))
    max_log2k = int(np.ceil(max(x)))
    ax.set_xticks(range(min_log2k, max_log2k + 1))
    ax.set_xticklabels(range(min_log2k, max_log2k + 1))

    ax.set_xlabel('log2(k)')
    
    y_min, y_max = ax.get_ylim()
    
    num_yticks = NUM_YTICKS  # You can adjust this number
    y_ticks = np.linspace(max(y_min, 1.0), y_max, num_yticks)
    # c_ticks = x_ticks * p
    y_labels = [f'{np.log2(y):.1f}' for y in y_ticks]

    plt.yticks(y_ticks, y_labels)
    
    ax.set_ylabel('Complexity (log2)')
    plt.title(f'Complexities for Succ. Prob. ≥ {success_probability} (m = 2{str(int(m_exp)).translate(translator)})')

    # transparent legend
    ax.legend(loc='upper left', fancybox=True, framealpha=0.5)
    # ax.legend()
    plt.grid(True)
    
    plt.subplots_adjust(bottom=0.2)

    plt.tight_layout()
    if pdf:
        pdf.savefig(bbox_inches='tight')
    elif output_path:
        plt.savefig(output_path, bbox_inches='tight')
    plt.close()

def main(args):
    csv_file = os.path.join(args.input_folder, 'results.csv')
    df = pd.read_csv(csv_file, comment='#')
    
    if args.output_mode == 'all':
        output_pdf = os.path.join(args.input_folder, 'plots.pdf')
        with PdfPages(output_pdf) as pdf:
            # Plot type 1
            for (m_exp, k), group in df.groupby(['m_exp', 'k']):
                if len(group) >= 0:
                    group = group.sort_values('n_coefficient')
                    plot_type1(group, m_exp, k, pdf)
                else:
                    print(f"Skipping group with m_exp={m_exp}, k={k} due to insufficient data points")
            
            # Plot type 2
            success_probabilities = [0.01, 0.5, 0.99]  # You can adjust these values as needed
            for m_exp, group in df.groupby('m_exp'):
                for success_probability in success_probabilities:
                    plot_type2(group, m_exp, success_probability, pdf)
    
    elif args.output_mode == 'ind':
        figures_folder = os.path.join(args.input_folder, 'figures')
        os.makedirs(figures_folder, exist_ok=True)
        
        # Plot type 1
        for (m_exp, k), group in df.groupby(['m_exp', 'k']):
            if len(group) >= 0:
                group = group.sort_values('n_coefficient')
                output_path = os.path.join(figures_folder, f'type1_m{int(m_exp)}_k{k}.pdf')
                plot_type1(group, m_exp, k, output_path=output_path)
            else:
                print(f"Skipping group with m_exp={m_exp}, k={k} due to insufficient data points")
        
        # Plot type 2
        success_probabilities = [0.01, 0.5, 0.99]  # You can adjust these values as needed
        for m_exp, group in df.groupby('m_exp'):
            for success_probability in success_probabilities:
                output_path = os.path.join(figures_folder, f'type2_m{int(m_exp)}_prob{success_probability}.pdf')
                plot_type2(group, m_exp, success_probability, output_path=output_path)
    
    elif args.output_mode == 'single':
        if args.plot_type == 'type1':
            group = df[(df['m_exp'] == args.m) & (df['k'] == args.k)]
            if len(group) > 0:
                output_path = os.path.join(args.input_folder, f'type1_m{args.m}_k{args.k}.pdf')
                plot_type1(group, args.m, args.k, output_path=output_path)
            else:
                print(f"No data found for m_exp={args.m} and k={args.k}")
                sys.exit(1)
        elif args.plot_type == 'type2':
            group = df[df['m_exp'] == args.m]
            if len(group) > 0:
                success_probabilities = [0.01, 0.5, 0.99]  # You can adjust these values as needed
                for success_probability in success_probabilities:
                    output_path = os.path.join(args.input_folder, f'type2_m{args.m}_prob{success_probability}.pdf')
                    plot_type2(group, args.m, success_probability, output_path=output_path)
            else:
                print(f"No data found for m_exp={args.m}")
                sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate plots from CSV data.")
    parser.add_argument("input_folder", help="Path to the input folder containing results.csv")
    parser.add_argument("-o", "--output-mode", choices=['all', 'ind', 'single'], required=True,
                        help="Output mode: 'all' for single PDF, 'ind' for individual PDFs, 'single' for a specific plot")
    parser.add_argument("-t", "--plot-type", choices=['type1', 'type2'],
                        help="Plot type for 'single' mode")
    parser.add_argument("-m", type=int, help="m_exp value for 'single' mode")
    parser.add_argument("-k", type=int, help="k value for 'single' mode with type1 plot")
    
    args = parser.parse_args()
    
    if args.output_mode == 'single':
        if not args.plot_type:
            parser.error("--plot-type is required when using -o single")
        if not args.m:
            parser.error("-m is required when using -o single")
        if args.plot_type == 'type1' and not args.k:
            parser.error("-k is required when using -o single with --plot-type type1")
    
    main(args)