import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, RepeatedStratifiedKFold
import os

BASE_DIR = "/Users/milkcaramelcheng/Desktop/project_1b"
MUTATION_FILE = os.path.join(BASE_DIR, "MFS_selected_genes", "MFS_selected_genes_mutation_Matrix.csv")
ANNOTATION_FILE = os.path.join(BASE_DIR, "MFS_selected_genes", "MFS_selected_genes.csv")
MIC_FILE = os.path.join(BASE_DIR, "MIC_data.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "MFS_selected_genes", "Elastic_Net_Outcomes_final")
OUTPUT_PLOT = os.path.join(OUTPUT_DIR, "ElasticNet_Feature_Importance_CHD.pdf")
OUTPUT_COMPARISON_PLOT = os.path.join(OUTPUT_DIR, "ElasticNet_L1_ratio_Comparison.pdf")
THRESHOLD_HIGH = 7.0
THRESHOLD_LOW = 4.0
MIN_MUTATION_COUNT = 3
L1_ratioS = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]

def parse_annotations(file_path):

    mapping = {}
    if not os.path.exists(file_path):
        print(f"Cannot find annotation file: {file_path}")
        return mapping
    try:
        df = pd.read_csv(file_path, header=None)
        annotation_col = 1 if df.shape[1] > 1 else None
        if annotation_col is None:
            print("Warning: No annotation column found, labels will be locus tags only.")
            return mapping
        for _, row in df.iterrows():
            tag = str(row[0]).strip()
            desc = str(row[annotation_col]).strip()
            if tag and tag != 'nan':
                if len(desc) > 50:
                    desc = desc[-50:].strip()
                mapping[tag] = desc
        print(f"Loaded {len(mapping)} annotations.")
    except Exception as e:
        print(f"Error reading annotation CSV: {e}")
    return mapping

def load_data_and_filter():
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
    X_raw = mut_df.loc[common]
    Y_raw = mic_df.loc[common]

    mutation_counts = X_raw.sum(axis=0)
    genes_to_keep = mutation_counts[mutation_counts >= MIN_MUTATION_COUNT].index
    X_filtered = X_raw[genes_to_keep]

    chd_col = next((c for c in Y_raw.columns if 'chd' in str(c).lower()), Y_raw.columns[0])
    mic_vals = Y_raw[chd_col]

    def categorize(val):
        v = float(val)
        if v >= THRESHOLD_HIGH: return "High"
        if v <= THRESHOLD_LOW: return "Low"
        return "Unknown"

    y_labels = mic_vals.apply(categorize)
    valid_idx = y_labels != "Unknown"

    return X_filtered.loc[valid_idx], y_labels[valid_idx]

def run_single_ratio(X_scaled, y_binary, genes, L1_ratio, C=0.5):

    model = LogisticRegression(
        l1_ratio=L1_ratio,
        C=C,
        solver='saga',
        class_weight='balanced',
        max_iter=5000,
        random_state=42
    )

    model.fit(X_scaled, y_binary)
    cv_strategy = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=42)
    cv_scores = cross_val_score(model, X_scaled, y_binary, cv=cv_strategy, scoring='accuracy')
    coefs = model.coef_[0]

    coef_df = pd.DataFrame({'Gene': genes, 'Coefficient': coefs})
    coef_df['Abs_Coefficient'] = coef_df['Coefficient'].abs()
    active = coef_df[coef_df['Abs_Coefficient'] > 1e-5].sort_values('Abs_Coefficient', ascending=False)

    return {
        'model': model,
        'cv_mean': cv_scores.mean(),
        'cv_std': cv_scores.std() / np.sqrt(len(cv_scores)),
        'n_active': len(active),
        'active_genes': active
    }

def plot_single_ratio(active_genes, L1_ratio, output_path):

    plt.rcParams['font.family'] = 'Arial'
    top_n = min(10, len(active_genes))
    if top_n == 0:
        print(f"L1 ratio={L1_ratio}: all coefficients shrunk to 0, skipping plot.")
        return
    plot_df = active_genes.head(top_n).copy()
    plot_df['Color'] = np.where(plot_df['Coefficient'] > 0, '#F73F43', '#60BFA4')

    fig_height = max(4, top_n * 0.35)
    plt.figure(figsize=(8, fig_height))
    sns.barplot(data=plot_df, x='Coefficient', y='Gene', hue='Gene',
                palette=plot_df['Color'].tolist(), width=0.7, legend=False)
    plt.axvline(x=0, color='black', linewidth=1)

    if L1_ratio == 0.0:
        ratio_label = "Pure Ridge (L1 ratio=0.0)"
    elif L1_ratio == 1.0:
        ratio_label = "Pure LASSO (L1 ratio=1.0)"
    else:
        ratio_label = f"Elastic Net (L1 ratio={L1_ratio})"

    plt.title(
        f"Top {top_n} Mutations Relevant to CHD Tolerance | {ratio_label}",
        fontsize=8, fontweight='bold'
    )
    plt.xlabel("Coefficient (Positive = Resistance, Negative = Susceptibility)", fontsize=6, fontweight='bold')
    plt.ylabel("", fontsize=6, fontweight='bold')
    plt.tick_params(axis='x', labelsize=6)
    plt.tick_params(axis='y', labelsize=6)
    ax = plt.gca()
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_comparison(summary_rows, output_path):

    plt.rcParams['font.family'] = 'Arial'
    df = pd.DataFrame(summary_rows)
    x_labels = [str(r) for r in df['L1_ratio']]

    fig, axes = plt.subplots(2, 1, figsize=(8, 8))

    # --- CV Accuracy ---
    axes[0].bar(x_labels, df['cv_mean'], yerr=df['cv_std'],
                color='#A0B2CF', capsize=3, width=0.7,
                error_kw=dict(elinewidth=0.8, ecolor='#333333', capthick=0.8))
    axes[0].set_ylim(0, 1.05)
    axes[0].set_title("Cross-Validation Accuracy vs L1 ratio", fontsize=12, fontweight='bold')
    axes[0].set_xlabel("L1 ratio  (0 = Ridge  →  1 = LASSO)", fontsize=8, fontweight='bold')
    axes[0].set_ylabel("CV Accuracy (mean ± SEM)", fontsize=8, fontweight='bold')
    axes[0].tick_params(labelsize=6)
    for label in axes[0].get_xticklabels() + axes[0].get_yticklabels():
        label.set_fontweight('bold')
    for i, (m, s) in enumerate(zip(df['cv_mean'], df['cv_std'])):
        axes[0].text(i, m + s + 0.01, f"{m:.2f}", ha='center', fontsize=7, fontweight='bold')
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)

    # --- Number of Active Features ---
    axes[1].bar(x_labels, df['n_active'], color='#EDAA6C', width=0.7)
    axes[1].set_title("Active Features (non-zero coef) vs L1 ratio", fontsize=12, fontweight='bold')
    axes[1].set_xlabel("L1 ratio  (0 = Ridge  →  1 = LASSO)", fontsize=8, fontweight='bold')
    axes[1].set_ylabel("Number of Active Genes", fontsize=8, fontweight='bold')
    axes[1].tick_params(labelsize=6)
    for label in axes[1].get_xticklabels() + axes[1].get_yticklabels():
        label.set_fontweight('bold')
    for i, n in enumerate(df['n_active']):
        axes[1].text(i, n + 0.3, str(n), ha='center', fontsize=7, fontweight='bold')
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)

    plt.suptitle("Elastic Net: Effect of LASSO vs Ridge Ratio on CHD Tolerance Model",
                 fontsize=11, y=1.02, fontweight='bold')
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.45)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def run_elastic_net_analysis():
    X, y_labels = load_data_and_filter()
    annotation_map = parse_annotations(ANNOTATION_FILE)

    y_binary = np.where(y_labels == 'High', 1, 0)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    genes = X.columns

    summary_rows = []
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for ratio in L1_ratioS:
        result = run_single_ratio(X_scaled, y_binary, genes, ratio)

        active = result['active_genes'].copy()
        active['Gene'] = active['Gene'].apply(
            lambda tag: f"{tag} ({annotation_map[tag]})" if tag in annotation_map else tag
        )
        result['active_genes'] = active

        penalty_label = "Ridge (L2)" if ratio == 0.0 else ("LASSO (L1)" if ratio == 1.0 else "ElasticNet")

        summary_rows.append({
            'L1_ratio': ratio,
            'penalty': penalty_label,
            'cv_mean': result['cv_mean'],
            'cv_std': result['cv_std'],
            'n_active': result['n_active']
        })

        ratio_str = str(ratio).replace('.', 'p')
        base, ext = os.path.splitext(OUTPUT_PLOT)
        ratio_plot_path = f"{base}_l1ratio_{ratio_str}{ext}"
        plot_single_ratio(result['active_genes'], ratio, ratio_plot_path)

    best = max(summary_rows, key=lambda r: r['cv_mean'])
    print(f"\n Best L1_ratio = {best['L1_ratio']} ({best['penalty']}) "
          f"with CV Accuracy = {best['cv_mean']:.2f}")

    plot_comparison(summary_rows, OUTPUT_COMPARISON_PLOT)

if __name__ == "__main__":
    run_elastic_net_analysis()