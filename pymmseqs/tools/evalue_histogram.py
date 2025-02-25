# pymmseqs/tools/evalue_histogram.py

import pandas as pd
import matplotlib.pyplot as plt
import os

from ..utils import has_header, to_superscript

def plot_evalues(
    input_file_path,
    start,
    end,
    step,
    output_dir,
    column_index=None,
):
    """
    Plots a histogram of E-values from MMseqs2 search results, automatically detecting the E-value column if possible.

    Args:
        input_file_path (str): Path to the tab-separated results file.
        column_index (int, optional): Index of the E-value column (0-based). If None, auto-detection is attempted.
        start (int): Starting exponent for E-value bins (e.g., -50 for 10⁻⁵⁰).
        end (int): Ending exponent for E-value bins (e.g., 0 for 10⁰).
        step (int): Step size between bin exponents (e.g., 5 for bins every 10⁵).
        output_dir (str): Directory to save the output plot.
    """
    # Read the file once with dtype=str to preserve original values
    df_str = pd.read_csv(input_file_path, sep='\t', header=None, dtype=str)

    # If column_index is provided, use it directly
    if column_index is not None:
        evalue_column = df_str.iloc[:, column_index]
    else:
        # Check if the file has a header using the provided has_header function
        if has_header(input_file_path):
            # Use the first row as column names
            header = df_str.iloc[0]
            df_str = df_str[1:]
            df_str.columns = header
            # Look for a column name indicating E-values
            possible_names = ['evalue', 'E-value', 'pvalue', 'P-value']
            evalue_col_name = None
            for col in df_str.columns:
                if str(col).lower() in possible_names:
                    evalue_col_name = col
                    break
            if evalue_col_name:
                evalue_column = df_str[evalue_col_name]
            else:
                # No matching column name; look for scientific notation
                evalue_column = None
                for col in df_str.columns:
                    if any('e-' in str(val).lower() for val in df_str[col].head(3)):
                        evalue_column = df_str[col]
                        break
                if evalue_column is None:
                    raise ValueError("Could not automatically detect the E-value column. Please specify column_index manually.")
        else:
            # No header; look for scientific notation in the first few rows
            evalue_column = None
            for j in range(len(df_str.columns)):
                if any('e-' in str(val).lower() for val in df_str.iloc[:, j].head(3)):
                    evalue_column = df_str.iloc[:, j]
                    break
            if evalue_column is None:
                raise ValueError("Could not automatically detect the E-value column. Please specify column_index manually.")

    # Convert E-value column to numeric, coercing invalid values to NaN
    evalue_column = pd.to_numeric(evalue_column, errors='coerce')

    # Define E-value thresholds
    thresholds = [10**i for i in range(start, end, step)]
    cumulative_counts = []
    non_cumulative_counts = []

    # Count sequences with E-values below each threshold
    for threshold in thresholds:
        cumulative_count = (evalue_column < threshold).sum()
        cumulative_counts.append(cumulative_count)

    # Calculate non-cumulative counts
    for i in range(len(cumulative_counts)):
        if i == 0:
            non_cumulative_counts.append(cumulative_counts[i])
        else:
            non_cumulative_counts.append(cumulative_counts[i] - cumulative_counts[i-1])

    # Create x-axis labels
    x_labels = list(range(start, end, step))
    bin_labels = []
    for i in range(len(x_labels)):
        if i == 0:
            bin_labels.append(f'≤10{to_superscript(x_labels[0])}')
        elif i == len(x_labels) - 1:
            bin_labels.append(f'≥10{to_superscript(x_labels[-1])}')
        else:
            bin_labels.append(f'10{to_superscript(x_labels[i])}')

    # Plotting
    plt.figure(figsize=(12, 9))
    plt.bar(range(len(bin_labels)), non_cumulative_counts, alpha=0.7, edgecolor='black', width=0.7)
    plt.xlabel('E-value threshold', fontsize=12)
    plt.ylabel('Number of sequences', fontsize=12)
    plt.title('Distribution of E-values in MMseqs2 Search Results', fontsize=14)
    plt.grid(axis='y', alpha=0.3, linestyle='--')
    plt.xticks(range(len(bin_labels)), bin_labels, rotation=45, ha='right', fontsize=11)
    plt.figtext(
        0.5, 0.01,
        f'Note: Bin 1 represents E-values ≤10{to_superscript(x_labels[0])}. Subsequent bins represent ranges from 10⁻ⁿ to 10⁻ᵐ, where n and m increase by {step}.\n'
        f'The last bin shows E-values ≥10{to_superscript(x_labels[-1])}. Lower E-values indicate higher significance.',
        ha='center', fontsize=10, style='italic'
    )
    plt.tight_layout(rect=[0, 0.05, 1, 0.97])

    # Save the plot
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "evalue_distribution.png")
    plt.savefig(output_path, dpi=300)
    print(f"Plot saved to: {output_path}")
    plt.close()
