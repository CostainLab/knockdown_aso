#### This script filters clinvar.vcf.gz downloaded from clinvar FTP 
#### for GrCH38 variants with sufficient evidence in GOF/DN gene list

# Load dependencies
import pandas as pd
import gzip
# import subprocess

# Input file
all_clinvar = '/hpf/projects/gregorycostain/kd_aso/clinvar/variant_summary/clinvar.vcf.gz'

# Output file
output = '/hpf/projects/gregorycostain/kd_aso/clinvar/variant_summary/clinvar_gof_12Jun26.tsv'

# Load gene list (mendeliome filtered for ASO-amenable) into a set
gene_list = '/hpf/projects/gregorycostain/kd_aso/metadata/gof_gene_list_g2p_logofunc.txt'

# Evidence level
review_status = ['practice_guideline', 'reviewed_by_expert_panel', 'criteria_provided,_multiple_submitters,_no_conflicts']

# Variant classification terms
pathogenic = r"\bPathogenic\b|\bLikely pathogenic\b"

with open(gene_list, "r") as f:
    # Removes newline characters from the end of each gene
    # Removes any blank lines
    # Ouputs as a set
    gene_set = {line.strip() for line in f if line.strip()}

# Helper: parse INFO column into dict
def parse_info(info_str):
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
    return info_dict


# Read VCF
rows = []
with gzip.open(all_clinvar, 'rt') as f:
    for line in f:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            header = line.strip().split('\t')
            continue
        
        values = line.strip().split('\t')
        row = dict(zip(header, values))
        
        info = parse_info(row['INFO'])
        
        # Extract relevant ClinVar fields
        gene = info.get('GENEINFO', '').split(':')[0]
        clin_sig = info.get('CLNSIG', '')
        review = info.get('CLNREVSTAT', '')
        origin = info.get('ORIGIN', '')
        variant_type = info.get('CLNVC', '')
        consequence = info.get('MC', '')
        
        rows.append({
            'CHROM': row['#CHROM'],
            'POS': row['POS'],
            'REF': row['REF'],
            'ALT': row['ALT'],
            'Type': variant_type,
            'GeneSymbol': gene,
            'ClinicalSignificance': clin_sig,
            'ReviewStatus': review,
            'Origin': origin,
            'Consequence': consequence
        })


# Convert to DataFrame
df = pd.DataFrame(rows)

# ----- Apply filters -----

# SNVs
df = df[df['Type'] == 'single_nucleotide_variant']

# Review status
df = df[df['ReviewStatus'].isin(review_status)]

# Gene
df = df[df['GeneSymbol'].isin(gene_set)]

# Germline
df = df[df['Origin'].str.contains(r'\b1\b', na=False, regex=True)]

# LPP
df = df[df['ClinicalSignificance'].str.contains(pathogenic, case=False, na=False)]

# Missense
df = df[df['Consequence'].str.contains('missense_variant', na=False)]


# Write output
df.to_csv(output, sep='\t', index=False)
# subprocess.run(['gzip', '-f', output], check=True)
