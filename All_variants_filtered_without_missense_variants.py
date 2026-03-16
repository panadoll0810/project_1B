import pandas as pd
import glob
import os

INPUT_FOLDER = "/Users/milkcaramelcheng/Desktop/project_1b/snippy/all_csv_files/*.csv"
OUTPUT_FOLDER = "/Users/milkcaramelcheng/Desktop/project_1b/All_variants_filtered_without_missense_variants"

def process_all_variants():

    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    snippy_files = glob.glob(INPUT_FOLDER)
    if not snippy_files:
        print(f"The Snippy CSV file was not found. Please check the path.")
        return

    for file_path in snippy_files:
        filename = os.path.basename(file_path)

        try:
            snp_df = pd.read_csv(file_path)
            if snp_df.empty:
                continue

            effect_col = 'EFFECT' if 'EFFECT' in snp_df.columns else 'effect'

            if effect_col in snp_df.columns:
                filtered_df = snp_df[~snp_df[effect_col].astype(str).str.contains('synonymous_variant|missense_variant', na=False)].copy()
            else:
                filtered_df = snp_df.copy()

            if filtered_df.empty:
                print(f"{filename}: No residual mutation after filtration.")
                continue

            output_path = os.path.join(OUTPUT_FOLDER, filename)
            filtered_df.to_csv(output_path, index=False)
            print(f"{filename}: saved {len(filtered_df)} mutations")

        except Exception as e:
            print(f"error while processing {filename} : {e}")

    print(f"done! outputs saved in {OUTPUT_FOLDER}")

if __name__ == "__main__":
    process_all_variants()