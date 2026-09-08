'''

Filtering functions for identifying putative histone modification regions.

'''

import numpy as np

import pandas as pd

ETA = 0.7

def load_chromosome_sizes(file_path):

    '''

    Load chromosome sizes from a two-column chromosome-size file.

    '''

    chromosome_sizes= pd.read_csv(file_path,sep='\t',header=None,names=['chromosome','chromosome_size'])

    return chromosome_sizes

def calculate_total_genomic_bins(chromosome_sizes,bin_size=1000):

    '''

    Calculate the total number of genomic bins across the chromosomes included in the analysis

    '''

    chromosome_sizes = chromosome_sizes.copy()

    chromosome_sizes['n_bins'] = np.ceil(chromosome_sizes['chromosome_size']/bin_size).astype(int)

    total_bins = chromosome_sizes['n_bins'].sum()

    return total_bins

def calculate_f_score(bin_counts,n_esc,n_npc):

    '''

    Calculate the combined normalized ChIP fragment enrichment score

    for each genomic bin

    '''

    bin_counts = bin_counts.copy()

    bin_counts['F_score'] = (bin_counts['esc_counts'] / n_esc) + (bin_counts['npc_counts'] / n_npc)

    return bin_counts

def calculate_putative_site_threshold(total_bins,eta=ETA):

    '''

    Calculate the ChIPDiff threshold used to identify putative histone modification sites

    '''

    threshold = 2/ (total_bins * eta)

    return threshold

def identify_putative_sites(bin_counts,threshold):

    '''

    Retain genomic bins whose combined ChIP signal exceeds the expected background threshold

    '''

    bin_counts = bin_counts.copy()

    bin_counts['putative_sites'] = (bin_counts['F_score'] > threshold)

    putative_sites = bin_counts[bin_counts['putative_sites']].copy().reset_index(drop=True)

    return putative_sites

def assign_region_ids(putative_sites):

    '''
    Group nearby putative modification sites into modification regions

    Consecutive putative bins on the same chromosome are assigned to the same region

    when the genomic gap between them at most 1kb
    '''

    putative_sites = putative_sites.sort_values(['chromosome','bin']).reset_index(drop=True)

    previous_chromosome = putative_sites['chromosome'].shift(1)

    previous_bin_end = putative_sites['bin_end'].shift(1)

    new_region = (putative_sites['chromosome'] != previous_chromosome) | (putative_sites['bin_end'] - previous_bin_end > 1000)

    putative_sites['region_id'] = new_region.cumsum().astype(int)

    return putative_sites

def build_putative_regions(bin_counts,chromosome_sizes,n_esc,n_npc):

    '''

    Calculate F-scores, identify putative modification sites and assign

    modification region IDs.

    '''

    total_bins = calculate_total_genomic_bins(chromosome_sizes)

    bin_counts = calculate_f_score(bin_counts,n_esc,n_npc)

    threshold = calculate_putative_site_threshold(total_bins)

    putative_sites = identify_putative_sites(bin_counts,threshold)

    putative_sites = assign_region_ids(putative_sites)

    return putative_sites