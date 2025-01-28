import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
from collections import Counter


def load_and_prepare_data(file_path):
    """Load data from a CSV file and prepare it for analysis"""
    # Load CSV file using comma as separator
    df = pd.read_csv(file_path, sep=',')

    # Drop the duplicate MotifLength column if it exists
    if 'FoundMotifLength' in df.columns and 'MotifLength' in df.columns:
        df = df.drop('FoundMotifLength', axis=1)

    # Parse MotifLengthDistribution into a more usable format
    if 'MotifLengthDistribution' in df.columns:
        df['MotifLengthDistribution'] = df['MotifLengthDistribution'].apply(
            lambda x: dict(pair.split(':') for pair in x.split())
        )

    return df


def plot_metrics_vs_threshold(df, motif_lengths):
    """Plot available metrics against threshold for each motif length"""
    # Define all possible metrics and check which ones are available
    all_metrics = {
        'GraphVertices': 'Number of Graph Vertices',
        'EdgeDensity': 'Edge Density',
        'ExecutionTimeMs': 'Execution Time (ms)',
        'MaxCliqueSize': 'Max Clique Size'
    }

    # Filter metrics to only those available in the DataFrame
    available_metrics = {k: v for k, v in all_metrics.items() if k in df.columns}

    if not available_metrics:
        raise ValueError("No valid metrics found in the DataFrame")

    # Calculate number of rows and columns for subplots
    n_metrics = len(available_metrics)
    n_cols = min(2, n_metrics)
    n_rows = (n_metrics + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 6 * n_rows))
    if n_metrics == 1:
        axes = [axes]
    else:
        axes = axes.ravel()

    colors = plt.cm.Set2(np.linspace(0, 1, len(motif_lengths)))

    for idx, (metric, title) in enumerate(available_metrics.items()):
        ax = axes[idx]
        for motif_len, color in zip(motif_lengths, colors):
            data = df[df['MotifLength'] == motif_len]
            ax.plot(data['Threshold'], data[metric],
                    marker='o', label=f'Length {motif_len}',
                    color=color, markersize=4)

        ax.set_xlabel('Threshold')
        ax.set_ylabel(title)
        ax.set_title(f'{title} vs Threshold')
        ax.grid(True, alpha=0.3)
        ax.legend()

    # Hide empty subplots if any
    for idx in range(len(available_metrics), len(axes)):
        axes[idx].set_visible(False)

    plt.tight_layout()
    return fig


def analyze_correlations(df):
    """Analyze correlations between metrics"""
    numeric_columns = ['Threshold', 'GraphVertices', 'EdgeDensity',
                       'ExecutionTimeMs', 'MaxCliqueSize']

    # Filter only the columns that exist in the DataFrame
    available_columns = [col for col in numeric_columns if col in df.columns]

    if not available_columns:
        print("No numeric columns available for correlation analysis")
        return None

    correlation_matrix = df[available_columns].corr()

    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0)
    plt.title('Correlation Matrix of Metrics')
    plt.tight_layout()

    return correlation_matrix


def analyze_threshold_effects(df):
    """Analyze how threshold affects various metrics"""
    results = {}
    metrics = ['GraphVertices', 'EdgeDensity', 'ExecutionTimeMs']

    # Filter only the metrics that exist in the DataFrame
    available_metrics = [metric for metric in metrics if metric in df.columns]

    for motif_len in df['MotifLength'].unique():
        motif_data = df[df['MotifLength'] == motif_len]
        correlations = {}

        for metric in available_metrics:
            correlation, p_value = stats.pearsonr(motif_data['Threshold'],
                                                  motif_data[metric])
            correlations[metric] = {
                'correlation': correlation,
                'p_value': p_value
            }

        results[motif_len] = correlations

    return results


def generate_summary_statistics(df):
    """Generate summary statistics for each motif length"""
    metrics = ['GraphVertices', 'EdgeDensity', 'ExecutionTimeMs', 'MaxCliqueSize']

    # Filter only the metrics that exist in the DataFrame
    available_metrics = [metric for metric in metrics if metric in df.columns]

    if not available_metrics:
        print("No metrics available for summary statistics")
        return None

    summary = df.groupby('MotifLength').agg({
        metric: ['mean', 'std', 'min', 'max'] for metric in available_metrics
    }).round(3)

    return summary


def plot_quality_parameters_heatmap(df):
    """Create a heatmap showing correlations between quality parameters"""
    quality_params = [
        'AverageMotifQuality', 'SequenceCoverage', 'TruePositiveRate',
        'MotifEntropy', 'AverageHammingDistance'
    ]

    # Filter only the columns that exist in the DataFrame
    available_params = [param for param in quality_params if param in df.columns]

    if not available_params:
        print("No quality parameters found in the data")
        return None

    correlation_matrix = df[available_params].corr()

    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, cmap='viridis', center=0)
    plt.title('Correlation Between Quality Parameters')
    plt.tight_layout()
    return plt.gcf()


def plot_motif_length_distribution(df):
    """Create a visualization of motif length distribution"""
    if 'MotifLengthDistribution' not in df.columns:
        print("Motif length distribution data not available")
        return None

    # Aggregate all distributions
    all_lengths = {}
    for dist_dict in df['MotifLengthDistribution']:
        for length, count in dist_dict.items():
            length = int(length)
            count = int(count)
            all_lengths[length] = all_lengths.get(length, 0) + count

    lengths = list(all_lengths.keys())
    counts = list(all_lengths.values())

    plt.figure(figsize=(10, 6))
    plt.bar(lengths, counts)
    plt.xlabel('Motif Length')
    plt.ylabel('Frequency')
    plt.title('Distribution of Found Motif Lengths')
    plt.grid(True, alpha=0.3)
    return plt.gcf()


def analyze_motif_similarity(df):
    """Visualize similarity between found motifs"""
    if 'Examples' not in df.columns:
        print("Motif examples not available in the data")
        return None

    def hamming_distance(s1, s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))

    # Extract unique motifs from the Examples column
    all_motifs = set()
    for example_str in df['Examples']:
        motifs = example_str.split('|')
        all_motifs.update(motifs)

    all_motifs = list(all_motifs)
    n_motifs = len(all_motifs)

    # Create similarity matrix
    similarity_matrix = np.zeros((n_motifs, n_motifs))
    for i, motif1 in enumerate(all_motifs):
        for j, motif2 in enumerate(all_motifs):
            if len(motif1) == len(motif2):
                similarity = 1 - hamming_distance(motif1, motif2) / len(motif1)
            else:
                similarity = 0
            similarity_matrix[i, j] = similarity

    plt.figure(figsize=(12, 10))
    sns.heatmap(similarity_matrix, xticklabels=all_motifs, yticklabels=all_motifs,
                cmap='YlOrRd', annot=True, fmt='.2f')
    plt.title('Motif Similarity Matrix')
    plt.tight_layout()
    return plt.gcf()


def plot_threshold_quality_relationship(df):
    """Visualize the relationship between threshold and motif quality metrics"""
    quality_metrics = ['AverageMotifQuality', 'SequenceCoverage',
                       'TruePositiveRate', 'MotifEntropy']

    # Filter only the metrics that exist in the DataFrame
    available_metrics = [metric for metric in quality_metrics if metric in df.columns]

    if not available_metrics:
        print("No quality metrics found in the data")
        return None

    n_metrics = len(available_metrics)
    fig, axes = plt.subplots(n_metrics, 1, figsize=(12, 4 * n_metrics))
    if n_metrics == 1:
        axes = [axes]

    for ax, metric in zip(axes, available_metrics):
        for motif_len in df['MotifLength'].unique():
            data = df[df['MotifLength'] == motif_len]
            ax.plot(data['Threshold'], data[metric],
                    marker='o', label=f'Length {motif_len}')

        ax.set_xlabel('Threshold')
        ax.set_ylabel(metric)
        ax.set_title(f'{metric} vs Threshold')
        ax.grid(True, alpha=0.3)
        ax.legend()

    plt.tight_layout()
    return fig


def analyze_motif_data(file_path):
    """Main analysis function for a given CSV file"""
    # Load and prepare data
    df = load_and_prepare_data(file_path)
    print(f"Analyzing file: {file_path}")
    print("Columns in DataFrame:", df.columns)

    try:
        # Original visualizations
        metrics_plot = plot_metrics_vs_threshold(df, sorted(df['MotifLength'].unique()))
        correlation_matrix = analyze_correlations(df)
        threshold_effects = analyze_threshold_effects(df)
        summary_stats = generate_summary_statistics(df)

        # New visualizations
        quality_heatmap = plot_quality_parameters_heatmap(df)
        length_dist = plot_motif_length_distribution(df)
        similarity_matrix = analyze_motif_similarity(df)
        threshold_quality = plot_threshold_quality_relationship(df)

        # Print analysis results
        print("\n=== Motif Analysis Results ===")
        print("\nSummary Statistics by Motif Length:")
        print(summary_stats)

        print("\nThreshold Effects Analysis:")
        for motif_len, results in threshold_effects.items():
            print(f"\nMotif Length {motif_len}:")
            for metric, stats_dict in results.items():
                print(f"  {metric}:")
                print(f"    Correlation with Threshold: {stats_dict['correlation']:.3f}")
                print(f"    P-value: {stats_dict['p_value']:.3f}")

        # Create distribution plots
        plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        df.boxplot(column='ExecutionTimeMs', by='MotifLength')
        plt.title('Execution Time Distribution by Motif Length')
        plt.suptitle('')

        plt.subplot(1, 2, 2)
        df.boxplot(column='GraphVertices', by='MotifLength')
        plt.title('Graph Vertices Distribution by Motif Length')
        plt.suptitle('')

        plt.tight_layout()
        plt.show()

    except Exception as e:
        print(f"Error during analysis: {str(e)}")


if __name__ == "__main__":
    target_directory = r"C:\Users\konrad_guest\Documents\Studia\semestr_3\Algorytmy_kombinatoryczne_w_bioinformatyce\zadanie_3\cmake-build-debug"

    print(f"Searching directory: {os.path.abspath(target_directory)}")

    if not os.path.isdir(target_directory):
        print("Directory does not exist:", target_directory)
    else:
        try:
            csv_files_pattern = os.path.join(target_directory, "*.csv")
            csv_files = glob.glob(csv_files_pattern)

            if not csv_files:
                print("No CSV files found in directory:", target_directory)
            else:
                print("Found CSV files:")
                for file in csv_files:
                    print(file)
                    analyze_motif_data(file)

        except Exception as e:
            print("Error reading directory:", e)