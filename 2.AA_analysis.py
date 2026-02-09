import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

BASE_DIR = "./data/integrated_analysis/"
INPUT_FILES = {
    "Cosmic": os.path.join(BASE_DIR, "Cosmic_with_ACDC_merged.tsv"),
    "Gnomad": os.path.join(BASE_DIR, "Gnomad_with_ACDC_merged.tsv")
}
OUTPUT_DIR = os.path.join(BASE_DIR, "transition_analysis")

AA_GROUP_MAP = {
    'G': 'Nonpolar', 'A': 'Nonpolar', 'V': 'Nonpolar', 'L': 'Nonpolar', 'I': 'Nonpolar',
    'M': 'Nonpolar', 'P': 'Nonpolar', 'F': 'Nonpolar', 'W': 'Nonpolar',
    'S': 'Polar', 'T': 'Polar', 'C': 'Polar', 'N': 'Polar', 'Q': 'Polar', 'Y': 'Polar',
    'K': 'Positive', 'R': 'Positive', 'H': 'Positive',
    'D': 'Negative', 'E': 'Negative'
}

CLASS_ORDER = ["Nonpolar", "Polar", "Negative", "Positive"]

def plot_heatmap(df, title, filename, label, cmap="Blues", is_float=False):
    """Generates a transition matrix heatmap."""
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(df.values, cmap=cmap)
    
    ax.set_title(title, fontsize=14, pad=15)
    ax.set_xticks(range(len(df.columns)))
    ax.set_xticklabels(df.columns, rotation=45, ha='right')
    ax.set_yticks(range(len(df.index)))
    ax.set_yticklabels(df.index)
    
    for (i, j), val in np.ndenumerate(df.values):
        text = f"{val:.1f}" if is_float else f"{int(val)}"
        ax.text(j, i, text, ha='center', va='center', fontsize=10)
    
    fig.colorbar(im, label=label)
    fig.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close(fig)

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    for name, path in INPUT_FILES.items():
        if not os.path.exists(path): continue
        print(f"Analyzing {name} transition matrix...")
        
        df = pd.read_csv(path, sep='\t', low_memory=False)
        
        if 'AA_ref_group' not in df.columns:
            df['AA_ref_group'] = df['AA_ref'].map(AA_GROUP_MAP)
            df['AA_alt_group'] = df['AA_alt'].map(AA_GROUP_MAP)
            df = df.dropna(subset=['AA_ref_group', 'AA_alt_group'])

        matrix = pd.crosstab(df['AA_ref_group'], df['AA_alt_group'])
        matrix = matrix.reindex(index=CLASS_ORDER, columns=CLASS_ORDER, fill_value=0)
        
        row_totals = matrix.sum(axis=1)
        percent_matrix = (matrix.T / row_totals).T * 100
        
        out_prefix = os.path.join(OUTPUT_DIR, f"{name.lower()}_transition")
        matrix.to_csv(f"{out_prefix}_counts.csv")
        percent_matrix.to_csv(f"{out_prefix}_percent.csv")
        
        plot_heatmap(matrix, f"{name} Variant Transition Counts", 
                     f"{out_prefix}_counts.png", "Count")
        plot_heatmap(percent_matrix, f"{name} Variant Transition (Row %)", 
                     f"{out_prefix}_percent.png", "Percentage", cmap="Greens", is_float=True)

    print("Transition analysis completed successfully.")

if __name__ == "__main__":
    main()