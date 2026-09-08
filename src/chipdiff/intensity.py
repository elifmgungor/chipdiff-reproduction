'''
Probabilistic histone-modification intensity estimation from ChipDIFF
'''
import numpy as np

ALPHA = 1
TAU = 3

def estimate_intensities(putative_sites,total_bins,n_esc, n_npc):
    '''
    Estimate posterior histone-modification intensities for ESC and NPC.
    A Beta(1,m) prior is used, where m is the total number of bins
    '''
    putative_sites = putative_sites.copy()
    beta_prior = total_bins
    putative_sites ['esc_intensity'] = (ALPHA + putative_sites['esc_counts']) / (ALPHA + beta_prior + n_esc)
    putative_sites ['npc_intensity'] = (ALPHA + putative_sites['npc_counts']) / (ALPHA + beta_prior +n_npc)
    return putative_sites

def calculate_fold_change(putative_sites, tau=TAU):
    '''
    Calculate ESC/NPC intensity ratios and provisional fold-change states.
    These states are descriptive and are not final HMM based ChIPDiff classifications.
    '''
    putative_sites = putative_sites.copy()
    putative_sites['intensity_ratio'] = putative_sites['esc_intensity'] / putative_sites['npc_intensity']
    putative_sites['log_intensity'] = np.log(putative_sites['intensity_ratio'])

    putative_sites['fold_change'] = 'non_differential'
    putative_sites.loc[putative_sites['intensity_ratio'] > tau, 'fold_change'] = 'ESC_enriched'
    putative_sites.loc[putative_sites['intensity_ratio'] < 1/tau, 'fold_change'] = 'NPC_enriched'
    return putative_sites

def calculate_posterior_parameters(putative_sites, total_bins, n_esc, n_npc):
    '''
    Calculate the Beta posterior parameters for ESC and NPC
    modification intensities in each bin
    '''
    putative_sites = putative_sites.copy()
    beta_prior = total_bins

    putative_sites['esc_alpha_post'] = ALPHA + putative_sites['esc_counts']
    putative_sites['esc_beta_post'] = beta_prior + n_esc - putative_sites['esc_counts']

    putative_sites['npc_alpha_post'] = ALPHA + putative_sites['npc_counts']
    putative_sites['npc_beta_post'] = beta_prior + n_npc - putative_sites['npc_counts']
    return putative_sites

def build_intensity_table(putative_sites, total_bins, n_esc, n_npc):
    '''
    Add posterior intensities, Beta posterior parameters, and provisional fold change
    classifications to putative sites
    '''
    putative_sites = estimate_intensities(putative_sites, total_bins, n_esc, n_npc)
    putative_sites = calculate_fold_change(putative_sites)
    putative_sites = calculate_posterior_parameters(putative_sites, total_bins, n_esc, n_npc)
    return putative_sites