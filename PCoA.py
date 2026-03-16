import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import os
from skbio.stats.distance import permanova, DistanceMatrix

BASE_DIR = "/Users/milkcaramelcheng/Desktop/project_1b"
MUTATION_FILE = os.path.join(BASE_DIR, "MFS_selected_genes", "MFS_regulator_mutation_Matrix.csv")
MIC_FILE = os.path.join(BASE_DIR, "MIC_data.csv")
OUTPUT_PLOT = os.path.join(BASE_DIR, "MFS_selected_genes", "PCoA_Jaccard_Clustering_MFS_regulator_with_label_test.pdf")
THRESHOLD_HIGH = 7.0
THRESHOLD_LOW = 4.0

def load_data():
    mut_df = pd.read_csv(MUTATION_FILE, index_col=0)
    strain_pattern = r'^(RS|PR|M\d|B\d)'
    cols_look_like_strains = mut_df.columns.astype(str).str.match(strain_pattern).sum()
    rows_look_like_strains = mut_df.index.astype(str).str.match(strain_pattern).sum()
    if cols_look_like_strains > rows_look_like_strains:
        mut_df = mut_df.T
    mut_df.index = mut_df.index.astype(str).str.replace(r'\.csv$', '', regex=True).str.strip()

    mic_df = pd.read_csv(MIC_FILE)
    strain_col = next((c for c in mic_df.columns if str(c).lower().strip() in ['strains', 'strain', 'sample']),
                      mic_df.columns[0])
    mic_df = mic_df.set_index(strain_col)
    mic_df.index = mic_df.index.astype(str).str.replace(r'\.csv$', '', regex=True).str.strip()

    common = mut_df.index.intersection(mic_df.index)
    if len(common) == 0: raise ValueError("Error - No common strains found between mutation matrix and MIC data. Please check the files and their index columns.")

    X = mut_df.loc[common]
    Y_df = mic_df.loc[common]

    chd_col = next((c for c in Y_df.columns if 'chd' in str(c).lower()), Y_df.columns[0])
    mic_vals = Y_df[chd_col]

    def categorize(val):
        try:
            v = float(val)
            if v >= THRESHOLD_HIGH:
                return "High"
            elif v <= THRESHOLD_LOW:
                return "Low"
            else:
                return "Intermediate"
        except:
            return "Unknown"

    groups = mic_vals.apply(categorize)
    return X, groups


def confidence_ellipse(x, y, ax, n_std=1.8, facecolor='none', **kwargs):
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    rx = np.sqrt(1 + pearson)
    ry = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=rx * 2, height=ry * 2,
                      facecolor=facecolor, **kwargs)
    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std
    transf = (transforms.Affine2D()
              .rotate_deg(45)
              .scale(scale_x, scale_y)
              .translate(np.mean(x), np.mean(y)))
    ellipse.set_transform(transf + ax.transData)
    ax.add_patch(ellipse)

def run_pcoa():

    try:
        X, groups = load_data()
    except Exception as e:
        print(e)
        return

    dist_matrix = squareform(pdist(X.values, metric='jaccard'))
    n = dist_matrix.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (dist_matrix ** 2) @ H
    eigvals, eigvecs = np.linalg.eigh(B)

    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]


    positive_eigvals = eigvals[eigvals > 0]
    total_variance = np.sum(positive_eigvals)
    explained_variance_ratio = positive_eigvals / total_variance


    pc1 = eigvecs[:, 0] * np.sqrt(eigvals[0])
    pc2 = eigvecs[:, 1] * np.sqrt(eigvals[1])

    pc1_var = explained_variance_ratio[0] * 100
    pc2_var = explained_variance_ratio[1] * 100

    dm = DistanceMatrix(dist_matrix, ids=list(X.index))
    grouping_series = pd.Series(groups.values, index=list(X.index), name='Group')
    permanova_result = permanova(dm, grouping_series, column='Group', permutations=999)

    np.random.seed(42)
    jitter_strength = 0.04 * (pc1.max() - pc1.min())
    pc1_jitter = pc1 + np.random.normal(0, jitter_strength, len(pc1))
    pc2_jitter = pc2 + np.random.normal(0, jitter_strength, len(pc2))

    plot_df = pd.DataFrame({
        'PCoA 1': pc1_jitter,
        'PCoA 2': pc2_jitter,
        'PCoA 1 orig': pc1,
        'PCoA 2 orig': pc2,
        'Group': groups.values,
        'Strain': X.index
    })

    fig, ax = plt.subplots(figsize=(11, 9))

    palette = {"High": "#F73F43", "Intermediate": "#ECBF51", "Low": "#60BFA4", "Unknown": "black"}
    markers = {"High": "o", "Intermediate": "o", "Low": "o", "Unknown": "o"}

    for group, color in palette.items():
        subset = plot_df[plot_df['Group'] == group]
        if len(subset) >= 3:
            confidence_ellipse(
                subset['PCoA 1'].values, subset['PCoA 2'].values, ax,
                n_std=1.8, edgecolor=color, linewidth=1.5,
                linestyle='--', alpha=0.6
            )

    for group in ["High", "Intermediate", "Low", "Unknown"]:
        subset = plot_df[plot_df['Group'] == group]
        if subset.empty:
            continue
        ax.scatter(
            subset['PCoA 1'], subset['PCoA 2'],
            label=group, color=palette[group],
            marker=markers[group],
            s=120, alpha=0.75, edgecolor='white', linewidth=0.8, zorder=3
        )

    from adjustText import adjust_text
    high_subset = plot_df[plot_df['Group'] == 'High']
    texts = []
    for _, row in high_subset.iterrows():
        texts.append(ax.text(
            row['PCoA 1'], row['PCoA 2'], row['Strain'],
            fontsize=7.5, color='#F73F43', fontweight='bold', zorder=5
        ))
    if texts:
        adjust_text(
            texts, ax=ax,
            autoalign=True,
            expand_points=(1.4, 1.4), expand_text=(1.2, 1.2),
            only_move={'points': 'xy', 'texts': 'xy'}
        )

    ax.set_xlabel(f"PCoA 1 ({pc1_var:.1f}%)", fontsize=12)
    ax.set_ylabel(f"PCoA 2 ({pc2_var:.1f}%)", fontsize=12)
    ax.set_title("PCoA of smvAR and tetAR (Jaccard Distance)", fontsize=14)

    p_val = permanova_result['p-value']
    pseudo_f = permanova_result['test statistic']
    p_str = f"p < 0.001" if p_val < 0.001 else f"p = {p_val:.3f}"
    ax.text(
        0.02, 0.02,
        f"PERMANOVA: pseudo-F = {pseudo_f:.3f}, {p_str}",
        transform=ax.transAxes, fontsize=9,
        verticalalignment='bottom', horizontalalignment='left',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7, edgecolor='gray')
    )

    ax.grid(False)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, title='CHD MIC Group', title_fontsize=10,
              fontsize=9, loc='best', framealpha=0.7)

    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300)
    print(f"All done! Saved in {OUTPUT_PLOT}")


if __name__ == "__main__":
    run_pcoa()