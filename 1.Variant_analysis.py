import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Path Configuration
BASE_DIR = "./data/integrated_analysis/"
GNOMAD_FILE = os.path.join(BASE_DIR, "Gnomad_with_ACDC_merged.tsv")
COSMIC_FILE = os.path.join(BASE_DIR, "Cosmic_with_ACDC_merged.tsv")

OUTPUT_DIR = os.path.join(BASE_DIR, "significant_variants")
OUT_STATS = os.path.join(OUTPUT_DIR, "ENSP_Specific_Thresholds.tsv")
OUT_SIG = os.path.join(OUTPUT_DIR, "Cosmic_ENSP_Specific_Significant.tsv")

# Threshold settings for outlier detection
LOWER_PCT = 2.5   
UPPER_PCT = 97.5  
MIN_VARIANTS = 10 

def main():
    """
    Identifies significant variants in COSMIC based on the null distribution
    derived from gnomAD for each specific protein (ENSP).
    """
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    print("Step 1: Loading gnomAD data for null distribution...")
    if not os.path.exists(GNOMAD_FILE):
        print(f"Error: {GNOMAD_FILE} not found.")
        return
        
    df_gnomad = pd.read_csv(GNOMAD_FILE, sep='\t', low_memory=False)
    
    if 'ENSP_core' not in df_gnomad.columns:
        if 'ENSP_or_Feature' in df_gnomad.columns:
            df_gnomad['ENSP_core'] = df_gnomad['ENSP_or_Feature'].astype(str).str.split('.').str[0]
        else:
            print("Error: ENSP information missing in gnomAD file.")
            return

    df_gnomad['ddg'] = pd.to_numeric(df_gnomad['ddg'], errors='coerce')
    df_gnomad = df_gnomad.dropna(subset=['ddg'])
    
    print("Step 2: Calculating thresholds per protein...")
    stats = df_gnomad.groupby('ENSP_core')['ddg'].agg(
        count='count',
        lower_bound=lambda x: np.percentile(x, LOWER_PCT),
        upper_bound=lambda x: np.percentile(x, UPPER_PCT)
    ).reset_index()
    
    stats['reliable'] = stats['count'] >= MIN_VARIANTS
    stats.to_csv(OUT_STATS, sep='\t', index=False)
    print(f" -> Analyzable proteins: {stats['reliable'].sum():,}")

    print("Step 3: Processing COSMIC data...")
    if not os.path.exists(COSMIC_FILE):
        print(f"Error: {COSMIC_FILE} not found.")
        return

    df_cosmic = pd.read_csv(COSMIC_FILE, sep='\t', low_memory=False)
    
    merged = pd.merge(df_cosmic, stats, on='ENSP_core', how='left')
    
    is_reliable = (merged['reliable'] == True)
    is_destabilizing = (merged['ddg'] < merged['lower_bound'])
    is_stabilizing = (merged['ddg'] > merged['upper_bound'])
    
    merged['Significance_Type'] = 'Not_Significant'
    merged.loc[is_reliable & is_destabilizing, 'Significance_Type'] = 'Extreme_Destabilizing'
    merged.loc[is_reliable & is_stabilizing, 'Significance_Type'] = 'Extreme_Stabilizing'
    
    significant_df = merged[merged['Significance_Type'] != 'Not_Significant'].copy()
    significant_df = significant_df.sort_values('ddg')
    significant_df.to_csv(OUT_SIG, sep='\t', index=False)
    
    print("\nSummary of Analysis:")
    print(f" - Total COSMIC variants: {len(df_cosmic):,}")
    print(f" - Significant variants: {len(significant_df):,} ({len(significant_df)/len(df_cosmic)*100:.2f}%)")
    print(f" - Results saved to: {OUT_SIG}")

if __name__ == "__main__":
    main()