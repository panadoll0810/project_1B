import pandas as pd
import glob
import os

BASE_DIR = "/Users/milkcaramelcheng/Desktop/project_1b"
INPUT_FOLDER = os.path.join(BASE_DIR, "All_variants_filtered_without_missense_variants", "*.csv")
GENE_LIST_FILE = os.path.join(BASE_DIR, "MFS_selected_genes", "MFS_selected_genes.csv")
OUTPUT_FILE = os.path.join(BASE_DIR, "MFS_selected_genes", "MFS_selected_genes_mutation_Matrix.csv")

def parse_gene_list(file_path):

    if not os.path.exists(file_path):
        print(f"The gene list file cannot be found: {file_path}")
        return []

    try:
        df = pd.read_csv(file_path, header=None, encoding='utf-8-sig')
        if df.empty:
            return []

        raw_values = df.iloc[:, 0].astype(str).str.strip().tolist()
        values = [v for v in raw_values if v and v.lower() not in {"nan", "none"}]

        if values:
            first = values[0].replace(" ", "").lower()
            if "locus" in first or "tag" in first:
                values = values[1:]

        tags = []
        seen = set()
        for tag in values:
            if tag not in seen:
                seen.add(tag)
                tags.append(tag)

        print(f"Loaded {len(tags)} locus tags from first column.")
        return tags
    except Exception as e:
        print(f"Error reading CSV gene list: {e}")
        return []

def generate_target_mutation_matrix():
    target_tags = parse_gene_list(GENE_LIST_FILE)

    if not target_tags:
        print("No gene IDs were extracted. Please check the content or path of the gene list file.")
        return

    files = glob.glob(INPUT_FOLDER)
    if not files:
        print("The CSV file of the variation results was not found. Please check the path.")
        return

    print(f"Found {len(files)} sample files，start screening and counting...")

    all_data = []
    all_filenames = []

    for file_path in files:
        filename = os.path.basename(file_path)
        all_filenames.append(filename)

        try:
            df = pd.read_csv(file_path)

            if df.empty:
                continue

            locus_col = 'LOCUS_TAG' if 'LOCUS_TAG' in df.columns else 'locus_tag'
            if locus_col not in df.columns:
                continue

            df[locus_col] = df[locus_col].astype(str)
            filtered_df = df[df[locus_col].isin(target_tags)]
            counts = filtered_df[locus_col].value_counts()

            for tag, count in counts.items():
                all_data.append({
                    'Locus_Tag': tag,
                    'Sample': filename,
                    'Count': count
                })

        except Exception as e:
            print(f"Error occurs while reading {filename} : {e}")

    if all_data:
        long_df = pd.DataFrame(all_data)
        matrix_df = long_df.pivot_table(index='Locus_Tag', columns='Sample', values='Count', fill_value=0)
    else:
        matrix_df = pd.DataFrame()

    matrix_df = matrix_df.reindex(columns=all_filenames, fill_value=0)
    matrix_df = matrix_df.reindex(index=target_tags, fill_value=0)
    matrix_df.columns = matrix_df.columns.str.replace('.csv', '', regex=False)
    matrix_df = matrix_df.sort_index(axis=1)
    matrix_df = matrix_df.sort_index(axis=0)

    matrix_df.to_csv(OUTPUT_FILE)
    print(f"All done! The mutation matrix for target genes has been generated.")
    print(f"Target gene (line): {len(matrix_df)}")
    print(f"Sample size (column): {len(matrix_df.columns)}")

if __name__ == "__main__":
    generate_target_mutation_matrix()