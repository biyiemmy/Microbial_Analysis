import pandas as pd
import os
import re

print("=== SIMPLE GTF ANNOTATION SCRIPT ===")

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

print(f"Processing genes: {gene_ids}")

# 4. Function to extract the best available info from GTF
def get_best_gtf_annotation(locus_tag):
    """Extract the best available annotation from GTF file"""
    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if locus_tag in line and not line.startswith('#'):
                    # Try to find any descriptive information
                    patterns = [
                        r'gene\s+"([^"]+)"',
                        r'product\s+"([^"]+)"',
                        r'note\s+"([^"]+)"',
                        r'description\s+"([^"]+)"',
                        r'db_xref\s+"([^"]+)"'
                    ]
                    
                    for pattern in patterns:
                        match = re.search(pattern, line)
                        if match:
                            return match.group(1)
                    
                    # If no patterns matched, return the whole attribute field
                    parts = line.split('\t')
                    if len(parts) > 8:
                        return parts[8][:150] + "..."  # Truncate if too long
                    
        return "No annotation found in GTF"
    except Exception as e:
        return f"Error: {str(e)}"

# 5. Get annotations
print("\n=== GENE ANNOTATIONS ===")
results = []

for gene_id in gene_ids:
    annotation = get_best_gtf_annotation(gene_id)
    print(f"{gene_id}: {annotation}")
    
    # Add statistical data
    gene_data = significant_degs[significant_degs.iloc[:, 0] == gene_id].iloc[0]
    results.append({
        "Locus_Tag": gene_id,
        "Annotation": annotation,
        "log2FoldChange": gene_data['logFC'],
        "Adjusted_P_value": gene_data['adj.P.Val'],
        "Average_Expression": gene_data['AveExpr']
    })

# 6. Save results
output_file = "results/Simple_Gene_Annotations.csv"
pd.DataFrame(results).to_csv(output_file, index=False)
print(f"\n✅ Simple annotations saved to {output_file}")

# 7. Print summary
print("\n" + "="*50)
print("TOP DIFFERENTIALLY EXPRESSED GENES")
print("="*50)

for result in results:
    print(f"\n{result['Locus_Tag']}")
    print(f"  FC: {result['log2FoldChange']:.2f}, p_adj: {result['Adjusted_P_value']:.2e}")
    print(f"  Annotation: {result['Annotation']}")