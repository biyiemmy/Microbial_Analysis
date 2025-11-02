import pandas as pd
import os
import re
import requests
import time

print("=== FINAL ENHANCED GTF + NCBI ANNOTATION SCRIPT ===")

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

# 4. Function to extract NCBI Gene ID from GTF
def get_ncbi_gene_id(locus_tag):
    """Extract NCBI Gene ID from GTF file"""
    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if locus_tag in line and not line.startswith('#'):
                    db_xref_match = re.search(r'db_xref\s+"GeneID:(\d+)"', line)
                    if db_xref_match:
                        return db_xref_match.group(1)
        return None
    except Exception as e:
        print(f"Error reading GTF for {locus_tag}: {e}")
        return None

# 5. Function to get gene info from NCBI
def get_ncbi_gene_info(ncbi_gene_id):
    """Get gene information from NCBI using Gene ID"""
    if not ncbi_gene_id:
        return "No NCBI Gene ID"
    
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={ncbi_gene_id}&retmode=json"
    
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            if 'result' in data and ncbi_gene_id in data['result']:
                gene_info = data['result'][ncbi_gene_id]
                name = gene_info.get('name', 'No name')
                description = gene_info.get('description', 'No description')
                summary = gene_info.get('summary', 'No summary')
                
                if summary != 'No summary':
                    return f"{name}: {summary[:200]}..."
                elif description != 'No description':
                    return f"{name}: {description}"
                else:
                    return f"{name} (no additional info)"
        
        return "Not found in NCBI"
    except Exception as e:
        return f"NCBI error: {str(e)}"

# 6. Get annotations
print("\n=== GENE ANNOTATIONS ===")
final_results = []

for i, gene_id in enumerate(gene_ids):
    print(f"Processing {gene_id} ({i+1}/{len(gene_ids)})...")
    
    # Get NCBI Gene ID from GTF
    ncbi_id = get_ncbi_gene_id(gene_id)
    print(f"  NCBI Gene ID: {ncbi_id}")
    
    # Get annotation from NCBI
    if ncbi_id:
        annotation = get_ncbi_gene_info(ncbi_id)
        print(f"  Annotation: {annotation}")
    else:
        annotation = "No NCBI Gene ID found in GTF"
        print(f"  {annotation}")
    
    # Add statistical data
    gene_data = significant_degs[significant_degs.iloc[:, 0] == gene_id].iloc[0]
    final_results.append({
        "Locus_Tag": gene_id,
        "NCBI_Gene_ID": ncbi_id,
        "Annotation": annotation,
        "log2FoldChange": gene_data['logFC'],
        "Adjusted_P_value": gene_data['adj.P.Val'],
        "Average_Expression": gene_data['AveExpr']
    })
    
    # Be polite to NCBI servers
    if i < len(gene_ids) - 1:
        time.sleep(1)

# 7. Save final results
output_file = "results/Final_Gene_Annotations.csv"
pd.DataFrame(final_results).to_csv(output_file, index=False)
print(f"\n✅ Final annotations saved to {output_file}")

# 8. Print beautiful summary
print("\n" + "="*60)
print("TOP DIFFERENTIALLY EXPRESSED GENES - FINAL ANNOTATIONS")
print("="*60)

for result in final_results:
    print(f"\n🔬 {result['Locus_Tag']}")
    print(f"   📊 Fold Change: {result['log2FoldChange']:.2f}")
    print(f"   📈 P-value: {result['Adjusted_P_value']:.2e}")
    print(f"   🧬 NCBI Gene ID: {result['NCBI_Gene_ID']}")
    print(f"   📝 Function: {result['Annotation']}")
    print("-" * 40)

print(f"\n🎉 Analysis complete! Results saved to: {output_file}")