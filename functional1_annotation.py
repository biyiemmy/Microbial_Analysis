import pandas as pd
import os
import re

# 1. Define paths
input_file = "results/DEGs_full_results.csv"
gtf_file = "ref_genome/Staphylococcus_aureus_ref.gtf"

# Check if files exist
if not os.path.exists(input_file):
    print(f"ERROR: {input_file} not found.")
    exit()
if not os.path.exists(gtf_file):
    print(f"ERROR: {gtf_file} not found.")
    exit()

# 2. Read DEG results
print("Reading DEG results...")
degs_df = pd.read_csv(input_file)
significant_degs = degs_df[degs_df['adj.P.Val'] < 0.05].copy()
significant_degs['abs_logFC'] = significant_degs['logFC'].abs()
significant_degs.sort_values('abs_logFC', ascending=False, inplace=True)

# 3. Get top genes
genes_to_annotate = significant_degs.head(10)
gene_ids = genes_to_annotate.iloc[:, 0].tolist()

print(f"Searching GTF for annotations: {gene_ids}")

# 4. Function to extract from GTF
def get_gtf_annotation(locus_tag):
    """Extract annotation from GTF file"""
    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                if locus_tag in line:
                    parts = line.split('\t')
                    if len(parts) > 8:
                        attributes = parts[8]
                        return f"Found: {attributes[:200]}..."  # Return first 200 chars
        return "Not found in GTF"
    except Exception as e:
        return f"Error: {str(e)}"

# 5. Get annotations
annotations = []
for gene_id in gene_ids:
    annotation = get_gtf_annotation(gene_id)
    annotations.append({"GeneID": gene_id, "Annotation": annotation})
    print(f"{gene_id}: {annotation}")

# 6. Save results
output_file = "results/GTF_Gene_Annotations.csv"
pd.DataFrame(annotations).to_csv(output_file, index=False)
print(f"\n✅ Saved GTF annotations to {output_file}")