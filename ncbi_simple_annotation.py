import pandas as pd
import os
import urllib.request
import time
from xml.etree import ElementTree as ET

print("=== NCBI GENE ANNOTATION SCRIPT (No requests needed) ===")

# 1. Define paths
input_file = "results/DEGs_full_results.csv"

if not os.path.exists(input_file):
    print(f"ERROR: {input_file} not found.")
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

# 4. Manual mapping based on NCBI Gene IDs we discovered
GENE_ID_MAPPING = {
    'SAOUHSC_T00048': '3921037',
    'SAOUHSC_T00034': '3921160', 
    'SAOUHSC_T00026': '3921152'
}

# 5. Function to get gene info from NCBI using urllib (no requests needed)
def get_ncbi_gene_info(ncbi_gene_id):
    """Get gene information from NCBI using Gene ID"""
    if not ncbi_gene_id:
        return "No NCBI Gene ID"
    
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={ncbi_gene_id}&retmode=xml"
    
    try:
        print(f"  Querying NCBI for Gene ID: {ncbi_gene_id}...")
        with urllib.request.urlopen(url) as response:
            xml_content = response.read().decode('utf-8')
            root = ET.fromstring(xml_content)
            
            # Try to extract gene description
            gene_desc = root.find('.//Gene-ref_desc')
            if gene_desc is not None and gene_desc.text:
                return gene_desc.text
            
            # Try to extract gene name
            gene_name = root.find('.//Gene-ref_locus')
            if gene_name is not None and gene_name.text:
                return gene_name.text
            
            return "No description available in XML response"
    except Exception as e:
        return f"Error: {str(e)}"

# 6. Get annotations from NCBI
print("\n=== GETTING ANNOTATIONS FROM NCBI ===")
final_results = []

for i, gene_id in enumerate(gene_ids):
    print(f"\nProcessing {gene_id} ({i+1}/{len(gene_ids)})...")
    
    # Get NCBI Gene ID from our mapping
    ncbi_id = GENE_ID_MAPPING.get(gene_id)
    
    if ncbi_id:
        print(f"  NCBI Gene ID: {ncbi_id}")
        annotation = get_ncbi_gene_info(ncbi_id)
        print(f"  Annotation: {annotation}")
    else:
        annotation = "No NCBI Gene ID mapping found"
        print(f"  {annotation}")
        ncbi_id = "N/A"
    
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
        time.sleep(2)

# 7. Save final results
output_file = "results/NCBI_Final_Annotations.csv"
pd.DataFrame(final_results).to_csv(output_file, index=False)
print(f"\n✅ Final NCBI annotations saved to {output_file}")

# 8. Print summary
print("\n" + "="*70)
print("TOP DIFFERENTIALLY EXPRESSED GENES - NCBI ANNOTATIONS")
print("="*70)

for result in final_results:
    print(f"\n{result['Locus_Tag']}")
    print(f"   Fold Change: {result['log2FoldChange']:.2f}")
    print(f"   P-value: {result['Adjusted_P_value']:.2e}")
    print(f"   NCBI Gene ID: {result['NCBI_Gene_ID']}")
    print(f"   Function: {result['Annotation']}")
    print("-" * 60)

print(f"\nAnalysis complete! Results saved to: {output_file}")