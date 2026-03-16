import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import math
import argparse

BASE_DIR = "/Users/milkcaramelcheng/Desktop/project_1b"
ORDER_FILE = os.path.join(BASE_DIR, "strain_order_CHD_MICs_hightolow.csv")
ANNOTATION_FILE = os.path.join(BASE_DIR, "MFS_selected_genes", "MFS_selected_genes.csv")
INPUT_MATRIX = os.path.join(BASE_DIR, "MFS_selected_genes", "MFS_selected_genes_mutation_Matrix.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "MFS_selected_genes", "Heatmap")

GENES_PER_PAGE = 100
CUSTOM_VMAX = 3
COLOR_MAP = mcolors.LinearSegmentedColormap.from_list("mic_heatmap",["#FAFBFD", "#A0B2CF", "#5A7AA8", "#2B3F5E"])
CELL_HEIGHT = 0.22
CELL_WIDTH = 0.22

def parse_csv_groups(file_path):

    if not os.path.exists(file_path):
        print(f"Cannot find file: {file_path}")
        return [], {}

    try:
        df = pd.read_csv(file_path)
        df.columns = [str(c).strip() for c in df.columns]

        for col in df.columns:
            df[col] = df[col].astype(str).str.strip()

        if df.empty:
            print("File is empty")
            return [], {}

        strain_col = 'Strains'
        if strain_col not in df.columns:
            strain_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]

        all_strains = df[strain_col].tolist()

        group_col = 'MIC range'

        if group_col in df.columns:
            df[group_col] = df[group_col].replace(['', 'nan', 'None'], None).ffill()
            strain_to_group = dict(zip(df[strain_col], df[group_col]))
        else:
            strain_to_group = {s: "Group1" for s in all_strains}

        return all_strains, strain_to_group

    except Exception as e:
        print(f"Failed to parse CSV: {e}")
        return [], {}


def parse_annotations(file_path):

    mapping = {}
    if not os.path.exists(file_path):
        print(f"Can not find annotation file: {file_path}")
        return mapping

    try:
        df = pd.read_csv(file_path, header=None)

        locus_col = 0
        annotation_col = 1 if df.shape[1] > 1 else None

        if annotation_col is None:
            print("Warning: No annotation column found, labels will be locus tags only.")
            return mapping

        for _, row in df.iterrows():
            tag = str(row[locus_col]).strip()
            desc = str(row[annotation_col]).strip()
            if tag and tag != 'nan':
                if len(desc) > 50:
                    desc = desc[-50:].strip()
                mapping[tag] = desc

        print(f"Loaded {len(mapping)-1} annotations from column index {annotation_col}.")
    except Exception as e:
        print(f"Error reading annotation CSV: {e}")

    return mapping

def plot_heatmap_final():
    if not os.path.exists(INPUT_MATRIX):
        print(f"Can not find the input matrix:{INPUT_MATRIX}")
        return

    ordered_samples, strain_group_map = parse_csv_groups(ORDER_FILE)
    if not ordered_samples:
        return

    annotation_map = parse_annotations(ANNOTATION_FILE)
    df = pd.read_csv(INPUT_MATRIX, index_col=0)
    df.columns = df.columns.astype(str).str.replace('.csv', '', regex=False).str.strip()
    available_samples = [s for s in ordered_samples if s in df.columns]

    if not available_samples:
        print("Error: No matching sample was found in the matrix.")
        return

    df_sorted = df.reindex(columns=ordered_samples, fill_value=0)
    df_sorted = df_sorted.fillna(0).astype(int)
    df_sorted = df_sorted.loc[(df_sorted != 0).any(axis=1)]

    new_index = []
    for tag in df_sorted.index:
        clean_tag = str(tag).strip()
        if clean_tag in annotation_map:
            new_label = f"{clean_tag} ({annotation_map[clean_tag]})"
        else:
            new_label = clean_tag
        new_index.append(new_label)
    df_sorted.index = new_index

    vlines = []
    group_spans = []
    if ordered_samples:
        current_g = strain_group_map.get(ordered_samples[0])
        seg_start = 0
        for i, sample in enumerate(ordered_samples):
            g = strain_group_map.get(sample)
            if g != current_g:
                vlines.append(i)
                group_spans.append((seg_start, i, current_g))
                seg_start = i
                current_g = g
        group_spans.append((seg_start, len(ordered_samples), current_g))

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    total_genes = len(df_sorted)
    total_pages = math.ceil(total_genes / GENES_PER_PAGE)

    for page in range(total_pages):
        start_idx = page * GENES_PER_PAGE
        end_idx = min((page + 1) * GENES_PER_PAGE, total_genes)
        subset_df = df_sorted.iloc[start_idx:end_idx]

        if subset_df.empty:
            continue

        fig_height = len(subset_df) * CELL_HEIGHT + 4
        fig_width = max(12, int(len(subset_df.columns) * CELL_WIDTH))

        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['font.size'] = 8
        plt.figure(figsize=(fig_width, fig_height))

        ax = sns.heatmap(
            subset_df,
            cmap=COLOR_MAP,
            linewidths=0.1,
            linecolor='white',
            annot=False,
            fmt="d",
            cbar_kws={'label': 'Mutation Count', 'shrink': 0.6},
            xticklabels=True,
            yticklabels=True,
            square=True,
            vmin=0,
            vmax=CUSTOM_VMAX
        )

        from matplotlib.ticker import MaxNLocator
        ax.collections[0].colorbar.locator = MaxNLocator(integer=True)
        ax.collections[0].colorbar.update_ticks()

        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_edgecolor('#333333')
            spine.set_linewidth(1.2)

        for x in vlines:
            ax.axvline(x, color='black', linewidth=2, linestyle='--')

        def get_mic_group_color(g_name):
            g_lower = str(g_name).lower()
            if 'high' in g_lower:
                return '#F73F43'
            elif 'low' in g_lower:
                return '#60BFA4'
            else:
                return '#ECBF51'

        group_colors = {g: get_mic_group_color(g) for _, _, g in group_spans}

        n_rows = len(subset_df)
        bar_height = max(0.5, n_rows * 0.03)
        bar_y_bottom = n_rows + 0.3

        for seg_start, seg_end, g_name in group_spans:
            color = group_colors.get(g_name, 'gray')

            ax.add_patch(plt.Rectangle(
                (seg_start, bar_y_bottom),
                seg_end - seg_start,
                bar_height,
                color=color,
                transform=ax.transData,
                clip_on=False,
                zorder=5
            ))

            ax.text(
                (seg_start + seg_end) / 2,
                bar_y_bottom + bar_height / 2,
                g_name,
                transform=ax.transData,
                ha='center', va='center',
                fontsize=8, fontweight='bold', color='white',
                clip_on=False,
                zorder=6
            )

        ax.xaxis.tick_top()
        plt.xticks(rotation=90, fontsize=8, fontweight='bold')
        plt.yticks(rotation=0, fontsize=8, fontweight='bold')
        plt.title(f"Mutations on MFS transporters and their putative regulators (Grouped by CHD MIC)", y=1.2, fontsize=14, fontweight='bold')
        out_name = f"Heatmap_Grouped_by_CHD_MIC_grouped.pdf"
        out_path = os.path.join(OUTPUT_DIR, out_name)
        plt.savefig(out_path, bbox_inches='tight', dpi=150)
        plt.close()

    print(f"All done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot paginated heatmaps from a mutation matrix.")
    parser.add_argument('--cell-height', type=float, default=CELL_HEIGHT,
                        help='Height of each heatmap cell in inches (default: 0.5)')
    parser.add_argument('--cell-width', type=float, default=CELL_WIDTH,
                        help='Width of each heatmap cell in inches (default: 0.25)')
    parser.add_argument('--genes-per-page', type=int, default=GENES_PER_PAGE,
                        help='Number of genes (rows) per page (default: 20)')
    args = parser.parse_args()

    CELL_HEIGHT = args.cell_height
    CELL_WIDTH = args.cell_width
    GENES_PER_PAGE = args.genes_per_page

    plot_heatmap_final()