'''
Post-processing functions for final ChIPDiff DHMS calls.

After HMM inference, genomic bins classified as ESC-enriched or
NPC-enriched are retained as differential histone modification sites (DHMSs).

Consecutive differential bins belonging to the same state are then 
merged into contunous genomic regions.
'''

import pandas as pd

DIFFERENTIAL_SITES = ['ESC_enriched','NPC_enriched']

def select_dhms_bins(hmm_table):
    '''
    Retain genomic bins classified as differential by the HMM
    A bin is considered differential when its final 'dhms_state' is
    either ESC-enriched or NPC-enriched.

    Parameters:
        -hmm_table: Bin-level DataFrame returned by HMM inference.
    Returns:
        -dhms_bins: DataFrame containing only differential bins. 
    '''
    
    dhms_bins = hmm_table[hmm_table['dhms_state'].isin(DIFFERENTIAL_SITES)].copy().reset_index(drop=True)
    return dhms_bins

def assign_dhms_region_ids(dhms_bins):
    '''
    Assign final region IDs to consecutive differential bins.
    A new DHMS region starts when:

    1. the chromosome changes 
    2. the differential state changes
    3. there is a genomic gap between consecutive 1kb bins.

    Therefore, only directly adjacent bins with the same differential state are merged.
    '''

    if dhms_bins.empty:
        dhms_bins = dhms_bins.copy()
        dhms_bins['dhms_region_id'] = pd.Series(dtype = int)
        return dhms_bins
    
    dhms_bins = dhms_bins.sort_values(["chromosome","bin_start"]).reset_index(drop=True)

    previous_chromosome = dhms_bins['chromosome'].shift(1)
    previous_bin_end = dhms_bins['bin_end'].shift(1)
    previous_state = dhms_bins['dhms_state'].shift(1)

    # Bins are considered continuous only when the current bin begins exactly 
    # where the previous bin ends. 

    new_region = ( (dhms_bins['chromosome'] != previous_chromosome)
                  | (dhms_bins['bin_start'] != previous_bin_end)
                  | (dhms_bins['dhms_state'] != previous_state)
                  )
    
    dhms_bins['dhms_region_id'] = new_region.cumsum().astype(int)

    return dhms_bins

def summarize_dhms_regions(dhms_bins):
    '''
    Collapse bin-level DHMS calls into continuous genomic regions.
    Returns one row per final DHMS region.
    '''
    if dhms_bins.empty:
        return pd.DataFrame(columns=["dhms_region_id","chromosome","region_start","region_end","dhms_state","n_bins",'mean_state_posterior'])
    
    # Store the posterior probability corresponding to the differential state assigned to each bin.
    dhms_bins = dhms_bins.copy()
    dhms_bins["state_posterior"] = dhms_bins['esc_posterior'].where(dhms_bins['dhms_state'] == "ESC_enriched", dhms_bins['npc_posterior'])

    dhms_regions = (dhms_bins.groupby("dhms_region_id",as_index=False)
                    .agg(chromosome=('chromosome','first'),region_start=('bin_start','min'),region_end=('bin_end','max'),
                    dhms_state=('dhms_state','first'),n_bins=('bin','size'),mean_state_posterior=('state_posterior','mean')) )
                
    return dhms_regions

def build_dhms_regions(hmm_table):
    """
    Build final ChIPDiff DHMS bin and region tables.
    Returns: 
    - dhms_bins : Differential genomic bins with final region IDs.
    - dhms_regions : Continuous DHMS regions obtained by merging adjacent
    differential bins of the same state.
    """

    dhms_bins = select_dhms_bins(hmm_table)
    dhms_bins = assign_dhms_region_ids(dhms_bins)
    dhms_regions = summarize_dhms_regions(dhms_bins)

    return (dhms_bins,dhms_regions)