import csv
import re
from Bio import SeqIO

GBK_FILE = "/Users/milkcaramelcheng/Desktop/project_1b/ref.gbk"
OUTPUT_CSV = "/Users/milkcaramelcheng/Desktop/project_1b/ref_gene_functions_with_protein_id.csv"

def extract_gene_info():
    try:
        records = SeqIO.parse(GBK_FILE, "genbank")
    except Exception as e:
        print(f"error: {e}")
        return

    with open(OUTPUT_CSV, "w", newline='', encoding='utf-8') as csvfile:
        fieldnames = ['Start', 'End', 'Strand', 'Gene_Symbol', 'Locus_Tag', 'Protein_ID', 'Product_Function']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()

        count = 0
        for record in records:
            for feature in record.features:
                if feature.type == "CDS":
                    start = int(feature.location.start) + 1
                    end = int(feature.location.end)
                    strand = "+" if feature.location.strand == 1 else "-"

                    gene = feature.qualifiers.get("gene", ["-"])[0]
                    locus = feature.qualifiers.get("locus_tag", ["-"])[0]
                    product = feature.qualifiers.get("product", ["Unknown Function"])[0]

                    raw_protein_id = feature.qualifiers.get("protein_id", ["-"])[0]
                    protein_match = re.search(r'(WP_\d+\.\d+)', raw_protein_id)
                    protein_id = protein_match.group(1) if protein_match else raw_protein_id

                    writer.writerow({
                        'Start': start,
                        'End': end,
                        'Strand': strand,
                        'Gene_Symbol': gene,
                        'Locus_Tag': locus,
                        'Protein_ID': protein_id,
                        'Product_Function': product
                    })
                    count += 1

    print(f"Done！Extracted {count} genes. Output saved to {OUTPUT_CSV}")


if __name__ == "__main__":
    extract_gene_info()