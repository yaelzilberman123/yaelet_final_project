import os
import random
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from codon_utils import sort_genes, dna_to_rna, count_codons, count_amino_acids, create_profile
from classify import classify_genes


def build_semi_profiles(semi_result_folder, codon_table):
    profiles = {}

    for file_name in sorted(os.listdir(semi_result_folder)):
        full_path = os.path.join(semi_result_folder, file_name)

        with open(full_path, 'r') as file:
            genes = sort_genes(file)

        random.shuffle(genes)
        split = int(len(genes) * 0.9)
        training_genes = genes[:split]
        test_genes = genes[split:]

        total_codon_counts = {codon: 0 for codon in codon_table}
        for gene in training_genes:
            rna = dna_to_rna(gene)
            gene_counts = count_codons(rna, codon_table)
            for codon in total_codon_counts:
                total_codon_counts[codon] += gene_counts[codon]

        amino_acid_counts = count_amino_acids(total_codon_counts, codon_table)
        profile = create_profile(total_codon_counts, amino_acid_counts, codon_table)
        profiles[file_name] = (profile, test_genes)

    return profiles


def generate_evo_bar_chart(accuracies, singular_name, output_folder):
    labels = [name.replace("_cds.fna", "") for name in accuracies]
    values = list(accuracies.values())

    plt.figure(figsize=(10, 5))
    bars = plt.bar(labels, values, color=[
        plt.cm.coolwarm(v / 100) for v in values
    ], edgecolor='black', linewidth=0.5)

    plt.ylim(0, 105)
    plt.ylabel("Classification Accuracy (%)")
    plt.xlabel("Comparison Organism")
    plt.title(f"Evolutionary Distance — {singular_name.replace('_cds.fna', '')}")
    plt.xticks(rotation=30, ha='right')
    plt.tight_layout()

    path = os.path.join(output_folder, singular_name.replace(".fna", "_bar.png"))
    plt.savefig(path, dpi=150)
    plt.show()


def generate_evo_heatmap(matrix, names, output_folder):
    labels = [n.replace("_cds.fna", "") for n in names]

    plt.figure(figsize=(10, 8))
    sns.heatmap(
        matrix,
        xticklabels=labels,
        yticklabels=labels,
        cmap="coolwarm",
        vmin=0, vmax=100,
        annot=True, fmt=".1f",
        linewidths=0.5
    )
    plt.title("Pairwise Classification Accuracy (%) — Evolutionary Distance")
    plt.xlabel("Comparison Organism")
    plt.ylabel("Singular Organism")
    plt.tight_layout()

    path = os.path.join(output_folder, "evo_dist_heatmap.png")
    plt.savefig(path, dpi=150)
    plt.show()

def generate_evo_scatter(accuracies, singular_name, evo_distances, output_folder, trendline_exclude=None):
    if trendline_exclude is None:
        trendline_exclude = []

    # drop self-comparison
    filtered = {k: v for k, v in accuracies.items() if k != singular_name}

    names = list(filtered.keys())
    x = [evo_distances[name] for name in names]
    y = [filtered[name] for name in names]
    labels = [name.replace("_cds.fna", "") for name in names]

    # trendline only on non-excluded points
    x_trend = [xi for xi, name in zip(x, names) if name not in trendline_exclude]
    y_trend = [yi for yi, name in zip(y, names) if name not in trendline_exclude]

    coeffs = np.polyfit(x_trend, y_trend, 1)
    trend = np.poly1d(coeffs)
    x_line = np.linspace(min(x_trend), max(x_trend), 200)

    y_pred = trend(np.array(x_trend))
    ss_res = np.sum((np.array(y_trend) - y_pred) ** 2)
    ss_tot = np.sum((np.array(y_trend) - np.mean(y_trend)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 1.0

    margin = (max(y) - min(y)) * 0.15
    y_bottom = min(y) - margin

    plt.figure(figsize=(10, 6))

    for xi, yi, label, name in zip(x, y, labels, names):
        if name in trendline_exclude:
            plt.scatter(xi, yi, color='tomato', s=80, zorder=3, marker='^')
            plt.annotate(f"{label}*", (xi, yi), textcoords="offset points", xytext=(8, 4), fontsize=9, color='tomato')
        else:
            plt.scatter(xi, yi, color='steelblue', s=80, zorder=3)
            plt.annotate(label, (xi, yi), textcoords="offset points", xytext=(8, 4), fontsize=9)

    plt.plot(x_line, trend(x_line), '--', color='gray', linewidth=1.5,
             label=f'Trendline (R² = {r_squared:.3f})')

    plt.xlabel("Evolutionary Distance (MYA)")
    plt.ylabel("Classification Accuracy (%)")
    plt.title(f"Codon Usage Accuracy vs Evolutionary Distance — {singular_name.replace('_cds.fna', '')}")
    plt.ylim(y_bottom, 105)
    plt.legend()

    if trendline_exclude:
        plt.figtext(0.99, 0.01, '* excluded from trendline (GC content outlier)',
                    ha='right', fontsize=8, color='tomato')

    plt.tight_layout()

    path = os.path.join(output_folder, singular_name.replace(".fna", "_scatter.png"))
    plt.savefig(path, dpi=150)
    plt.show()

def run_evo_dist(
    singular_file_name,
    codon_table,
    evo_distances,
    trendline_exclude=None,
    semi_result_folder="results/semi_result",
    output_folder="results/evodist_results"
):
    profiles = build_semi_profiles(semi_result_folder, codon_table)
    singular_profile, test_genes = profiles[singular_file_name]

    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, singular_file_name.replace(".fna", "_evo_dist.txt"))

    accuracies = {}

    with open(output_path, 'w') as out_file:
        header = f"---- Evolutionary distance classification: {singular_file_name} ----\n"
        print(header, end="")
        out_file.write(header)

        for other_name, (other_profile, _) in sorted(profiles.items()):
            accuracy = classify_genes(
                test_genes,
                [singular_profile, other_profile],
                singular_profile,
                codon_table
            )
            accuracies[other_name] = accuracy
            line = f"{other_name}: {accuracy:.2f}%\n"
            print(line, end="")
            out_file.write(line)

    generate_evo_scatter(accuracies, singular_file_name, evo_distances, output_folder, trendline_exclude)

def run_full_evo_dist(
    codon_table,
    semi_result_folder="results/semi_result",
    output_folder="results/evodist_results"
):
    profiles = build_semi_profiles(semi_result_folder, codon_table)
    names = sorted(profiles.keys())
    n = len(names)
    matrix = np.zeros((n, n))

    os.makedirs(output_folder, exist_ok=True)

    for i, singular_name in enumerate(names):
        singular_profile, test_genes = profiles[singular_name]
        for j, other_name in enumerate(names):
            other_profile, _ = profiles[other_name]
            matrix[i][j] = classify_genes(
                test_genes,
                [singular_profile, other_profile],
                singular_profile,
                codon_table
            )

    generate_evo_heatmap(matrix, names, output_folder)