#!/usr/bin/env python3
"""
Immune Gene Reference Builder
Downloads and formats immune gene coordinates for Fst analysis
"""

import pandas as pd
import requests
import json
from io import StringIO
import time

def download_ensembl_genes(gene_symbols, build='GRCh38'):
    """
    Download gene coordinates from Ensembl REST API
    """
    base_url = "https://rest.ensembl.org"
    
    results = []
    
    for gene in gene_symbols:
        try:
            # Get gene info
            url = f"{base_url}/lookup/symbol/homo_sapiens/{gene}"
            headers = {"Content-Type": "application/json"}
            
            response = requests.get(url, headers=headers)
            
            if response.status_code == 200:
                data = response.json()
                
                results.append({
                    'gene_symbol': gene,
                    'ensembl_id': data.get('id', ''),
                    'chrom': data.get('seq_region_name', ''),
                    'start': data.get('start', 0),
                    'end': data.get('end', 0),
                    'strand': data.get('strand', 0),
                    'biotype': data.get('biotype', ''),
                    'description': data.get('description', '')
                })
                
                print(f"Downloaded: {gene}")
                
            else:
                print(f"Failed to download: {gene}")
                
            # Rate limiting
            time.sleep(0.1)
            
        except Exception as e:
            print(f"Error downloading {gene}: {e}")
    
    return pd.DataFrame(results)

def get_immune_gene_sets():
    """
    Define comprehensive immune gene sets
    """
    
    # HLA genes
    hla_genes = [
        'HLA-A', 'HLA-B', 'HLA-C',  # Class I
        'HLA-DRA', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5',  # Class II DR
        'HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'HLA-DQB2',  # Class II DQ
        'HLA-DPA1', 'HLA-DPB1', 'HLA-DPB2',  # Class II DP
        'HLA-E', 'HLA-F', 'HLA-G',  # Non-classical Class I
        'HLA-DOA', 'HLA-DOB', 'HLA-DMA', 'HLA-DMB'  # Class II processing
    ]
    
    # Toll-like receptors
    tlr_genes = [
        'TLR1', 'TLR2', 'TLR3', 'TLR4', 'TLR5', 'TLR6',
        'TLR7', 'TLR8', 'TLR9', 'TLR10'
    ]
    
    # Cytokines and chemokines
    cytokine_genes = [
        'IL1A', 'IL1B', 'IL1RN', 'IL2', 'IL4', 'IL5', 'IL6', 'IL8',
        'IL10', 'IL12A', 'IL12B', 'IL13', 'IL15', 'IL17A', 'IL18',
        'TNF', 'TNFA', 'IFNG', 'IFNA1', 'IFNB1',
        'CCL2', 'CCL3', 'CCL4', 'CCL5', 'CXCL8', 'CXCL10',
        'CCR1', 'CCR2', 'CCR3', 'CCR4', 'CCR5', 'CCR7',
        'CXCR1', 'CXCR2', 'CXCR3', 'CXCR4', 'CXCR5'
    ]
    
    # Complement system
    complement_genes = [
        'C1QA', 'C1QB', 'C1QC', 'C1R', 'C1S',
        'C2', 'C3', 'C4A', 'C4B', 'C5', 'C6', 'C7', 'C8A', 'C8B', 'C8G', 'C9',
        'CFB', 'CFD', 'CFH', 'CFI', 'CFP',
        'CR1', 'CR2', 'CR3', 'CD55', 'CD46', 'CD59'
    ]
    
    # Immunoglobulins
    ig_genes = [
        'IGHV1-69', 'IGHV3-23', 'IGHV4-34',  # Heavy chain variable
        'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4',  # Heavy chain constant
        'IGHA1', 'IGHA2', 'IGHM', 'IGHD', 'IGHE',
        'IGKV1-39', 'IGKV2-28', 'IGKV3-20',  # Kappa light chain
        'IGLV1-44', 'IGLV2-14', 'IGLV3-25'   # Lambda light chain
    ]
    
    # Antimicrobial peptides
    amp_genes = [
        'DEFB1', 'DEFB4A', 'DEFB103A', 'DEFB104A', 'DEFB105A',
        'DEFA1', 'DEFA3', 'DEFA4', 'DEFA5', 'DEFA6',
        'CAMP', 'LCN2', 'LYZ', 'SLPI'
    ]
    
    # Antiviral genes
    antiviral_genes = [
        'MX1', 'MX2', 'OAS1', 'OAS2', 'OAS3', 'OASL',
        'IFITM1', 'IFITM2', 'IFITM3', 'IFITM5',
        'PKR', 'RIG1', 'MDA5', 'MAVS', 'TRIM25'
    ]
    
    # Pathogen-specific resistance genes
    pathogen_resistance = [
        'ACKR1',  # Duffy (malaria)
        'G6PD',   # Glucose-6-phosphate dehydrogenase (malaria)
        'HBB',    # Beta-globin (sickle cell, malaria)
        'HBA1', 'HBA2',  # Alpha-globin (thalassemia, malaria)
        'SLC11A1',  # NRAMP1 (tuberculosis, leishmaniasis)
        'VDR',    # Vitamin D receptor (tuberculosis)
        'ERAP2',  # Endoplasmic reticulum aminopeptidase 2 (plague)
        'FUT2',   # Fucosyltransferase 2 (norovirus, cholera)
        'APOL1'   # Apolipoprotein L1 (trypanosomiasis)
    ]
    
    # Combine all gene sets
    all_immune_genes = (hla_genes + tlr_genes + cytokine_genes + 
                       complement_genes + ig_genes + amp_genes + 
                       antiviral_genes + pathogen_resistance)
    
    return {
        'HLA': hla_genes,
        'TLR': tlr_genes,
        'Cytokines': cytokine_genes,
        'Complement': complement_genes,
        'Immunoglobulins': ig_genes,
        'Antimicrobial_Peptides': amp_genes,
        'Antiviral': antiviral_genes,
        'Pathogen_Resistance': pathogen_resistance,
        'All_Immune': list(set(all_immune_genes))  # Remove duplicates
    }

def create_immune_gene_reference(output_file='immune_genes_reference.csv', 
                                build='GRCh38'):
    """
    Create comprehensive immune gene reference file
    """
    
    print("Getting immune gene sets...")
    gene_sets = get_immune_gene_sets()
    
    print(f"Downloading coordinates for {len(gene_sets['All_Immune'])} genes...")
    
    # Download gene coordinates
    immune_genes_df = download_ensembl_genes(gene_sets['All_Immune'], build)
    
    # Add gene categories
    immune_genes_df['category'] = ''
    
    for category, genes in gene_sets.items():
        if category != 'All_Immune':
            mask = immune_genes_df['gene_symbol'].isin(genes)
            immune_genes_df.loc[mask, 'category'] = (
                immune_genes_df.loc[mask, 'category'] + 
                category + ';'
            )
    
    # Clean up category column
    immune_genes_df['category'] = immune_genes_df['category'].str.rstrip(';')
    
    # Filter valid chromosomes (remove patches, alternative sequences)
    valid_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
    immune_genes_df = immune_genes_df[
        immune_genes_df['chrom'].isin(valid_chroms)
    ]
    
    # Sort by chromosome and position
    immune_genes_df['chrom_num'] = immune_genes_df['chrom'].replace({
        'X': 23, 'Y': 24, 'MT': 25
    })
    immune_genes_df['chrom_num'] = pd.to_numeric(
        immune_genes_df['chrom_num'], errors='coerce'
    ).fillna(26)
    
    immune_genes_df = immune_genes_df.sort_values(['chrom_num', 'start'])
    immune_genes_df = immune_genes_df.drop('chrom_num', axis=1)
    
    # Save to file
    immune_genes_df.to_csv(output_file, index=False)
    
    print(f"Created immune gene reference: {output_file}")
    print(f"Total genes: {len(immune_genes_df)}")
    print(f"Chromosomes: {sorted(immune_genes_df['chrom'].unique())}")
    
    # Print summary by category
    print("\nGenes by category:")
    for category, genes in gene_sets.items():
        if category != 'All_Immune':
            count = len(set(genes) & set(immune_genes_df['gene_symbol']))
            print(f"  {category}: {count}/{len(genes)} genes found")
    
    return immune_genes_df

def download_immport_genelist(list_name='innate_immune_genes'):
    """
    Download gene list from ImmPort (requires manual download)
    This function shows how to format ImmPort data
    """
    
    # ImmPort gene lists need to be downloaded manually from:
    # https://www.immport.org/shared/geneList
    
    # Example of formatting downloaded ImmPort data
    immport_example = """
    # After downloading from ImmPort, format like this:
    
    Gene Symbol,Gene Name,Category
    TLR1,toll like receptor 1,Innate
    TLR2,toll like receptor 2,Innate
    TLR3,toll like receptor 3,Innate
    """
    
    print("To download ImmPort gene lists:")
    print("1. Go to https://www.immport.org/shared/geneList")
    print("2. Select desired gene list (e.g., 'Innate immune genes')")
    print("3. Download as CSV")
    print("4. Use the format_immport_data() function to process")
    
    return immport_example

def format_immport_data(immport_file):
    """
    Format downloaded ImmPort gene list
    """
    
    # Read ImmPort CSV
    immport_df = pd.read_csv(immport_file)
    
    # Standardize column names
    immport_df.columns = immport_df.columns.str.lower().str.replace(' ', '_')
    
    # Get coordinates from Ensembl
    gene_symbols = immport_df['gene_symbol'].tolist()
    coordinates_df = download_ensembl_genes(gene_symbols)
    
    # Merge with ImmPort data
    merged_df = pd.merge(
        coordinates_df, 
        immport_df, 
        left_on='gene_symbol', 
        right_on='gene_symbol',
        how='left'
    )
    
    return merged_df

# Roman period specific immune genes of interest
def get_roman_period_immune_genes():
    """
    Immune genes particularly relevant for Roman period analysis
    """
    
    roman_immune_genes = {
        'Urban_Pathogens': [
            'TLR2', 'TLR4', 'TLR5',  # Bacterial recognition
            'IL1B', 'IL6', 'TNF',   # Inflammatory response
            'CCR5', 'CCR2',         # Chemokine receptors
            'DEFB1', 'DEFB4A'       # Antimicrobial peptides
        ],
        
        'Malaria_Resistance': [
            'ACKR1', 'G6PD', 'HBB', 'HBA1', 'HBA2'
        ],
        
        'Tuberculosis_Resistance': [
            'SLC11A1', 'VDR', 'IL12B', 'IFNG'
        ],
        
        'Plague_Resistance': [
            'CCR5', 'ERAP2', 'TLR4'
        ],
        
        'Dietary_Adaptation': [
            'FUT2',   # Gut microbiome
            'VDR',    # Vitamin D (latitude adaptation)
            'APOA1', 'APOB'  # Lipid metabolism
        ],
        
        'HLA_Diversity': [
            'HLA-A', 'HLA-B', 'HLA-C',
            'HLA-DRB1', 'HLA-DQB1', 'HLA-DPB1'
        ]
    }
    
    return roman_immune_genes

if __name__ == "__main__":
    # Create comprehensive immune gene reference
    immune_genes_df = create_immune_gene_reference()
    
    # Also create Roman period specific reference
    roman_genes = get_roman_period_immune_genes()
    all_roman_genes = []
    for category, genes in roman_genes.items():
        all_roman_genes.extend(genes)
    
    roman_genes_df = download_ensembl_genes(list(set(all_roman_genes)))
    roman_genes_df.to_csv('roman_period_immune_genes.csv', index=False)
    
    print("\nCreated files:")
    print("- immune_genes_reference.csv (comprehensive)")
    print("- roman_period_immune_genes.csv (Roman period specific)")
