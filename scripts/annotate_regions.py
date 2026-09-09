from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
REGIONS_FILE = PROJECT_ROOT/"results/tables/dhms_regions.csv"
REFSEQ_FILE = PROJECT_ROOT/"metadata/mm8.refFlat.txt.gz"
OUTPUT_FILE = PROJECT_ROOT/"results/tables/dhms_regions_annotated.csv"
PROMOTER_DISTANCE = 1000

REFLAT_COLUMNS = ["gene_name","transcript_id","chromosome","strand","tx_start",
                 "tx_end","cds_start","cds_end","exon_count","exon_starts","exon_ends"]


def load_refseq_genes():
    '''
    Load UCSC RefSeq annotations for the mm8 mouse genome assembly.

    When multiple transcripts map to the same gene, retain the transcript
    with the longest coding sequence, following the strategy described in the original
    ChIPDiff study
    '''

    genes = pd.read_csv(REFSEQ_FILE,names=REFLAT_COLUMNS, sep='\t')
    canonical_chromosomes = [f"chr{id}" for id in range(1,20)] + ['chrX', "chrY"]
    genes = genes[genes['chromosome'].isin(canonical_chromosomes)].copy()

    #Coding sequence lengths are used to select one transcript per gene.
    genes['orf_length'] = genes['cds_end'] -genes['cds_start']  
    genes = genes.sort_values('orf_length',ascending=False).drop_duplicates(
            subset="gene_name",keep="first").reset_index(drop=True)

    #Transcription start site depends on gene strand
    genes['tss'] = np.where(genes["strand"]== "+", genes["tx_start"],genes["tx_end"])
    #Promoters are defined as +/- 1kb around the TSS
    genes['promoter_start'] = (genes['tss']-PROMOTER_DISTANCE).clip(lower=0)
    genes['promoter_end'] = genes['tss'] + PROMOTER_DISTANCE
    
    return genes

def annotate_regions(regions,genes):
    '''
    Annotate DHMS regions chromosome by chromosome
    '''
    annotated_regions = []

    for chromosome,chromosome_regions in regions.groupby("chromosome",sort=False):
        #Keep only RefSeq genes located on the chromosome currently being annotated
        chromosome_genes = genes[genes['chromosome']==chromosome]
        annotations = chromosome_regions.apply(annotate_region, axis=1,chromosome_genes=chromosome_genes)
        chromosome_result = pd.concat([chromosome_regions.reset_index(drop=True),annotations.reset_index(drop=True)],axis=1)
        annotated_regions.append(chromosome_result)
    return pd.concat(annotated_regions,ignore_index=True)


def annotate_region(region,chromosome_genes):
    '''
    Annotate one DHMS region relative to RefSeq genes.
    Annotation priority: promoter > gene_body > intergenic
    The nearest gene and distance to its TSS are also reported
    '''

    region_start = region['region_start']
    region_end = region['region_end']
    region_center = (region_start + region_end) / 2

    # First check: does DHMS overlap a promoter
    promoter_overlap = chromosome_genes[(chromosome_genes['promoter_start'] <region_end) & (chromosome_genes['promoter_end']>region_start)]

    if not promoter_overlap.empty:
        distances = (promoter_overlap['tss'] - region_center).abs()
        # idxmin finds the index of the smallest distance
        # loc. uses this index to take all the info
        gene= promoter_overlap.loc[distances.idxmin()]
        annotation = "promoter"
    
    else: 
        # If no promoter overlaps, test gene body overlap
        gene_body_overlap = chromosome_genes[(chromosome_genes['tx_start']<region_end) & (chromosome_genes['tx_end'] > region_start)]
        if not gene_body_overlap.empty:
            distances = (gene_body_overlap['tss']- region_center).abs()
            gene = gene_body_overlap.loc[distances.idxmin()]
            annotation = "gene_body"
        else:
            gene=None
            annotation = "intergenic"
    
    #Find the closest gene independently of overlap
    nearest_index = (chromosome_genes['tss'] - region_center).abs().idxmin()
    nearest_gene = chromosome_genes.loc[nearest_index]
    distance_to_tss = region_center - nearest_gene['tss']

    if gene is None:
        gene_name = np.nan
        transcript_id = np.nan
        strand = np.nan
    else:
        gene_name = gene['gene_name']
        transcript_id = gene['transcript_id']
        strand = gene['strand']
    
    return pd.Series({"annotation":annotation,"gene_name":gene_name,"transcript_id":transcript_id,
                    "strand":strand,"nearest_gene":nearest_gene['gene_name'],"distance_to_tss":distance_to_tss})
    
def main():
    print("Loading DHMS regions...")
    regions = pd.read_csv(REGIONS_FILE)
    print("Loading mm8 RefSeq genes...")
    genes = load_refseq_genes()
    print(f"Unique RefSeq genes retained {len(genes)}")
    annotated_regions = annotate_regions(regions,genes)
    annotated_regions.to_csv(OUTPUT_FILE,index=False)

    print("Genomic annotation summary")
    print("--------------------------")
    print(annotated_regions['annotation'].value_counts())
    print("Annotations by HMM state")
    print(pd.crosstab(annotated_regions["annotation"],annotated_regions["dhms_state"]))
    print(f"Annotated regions saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()