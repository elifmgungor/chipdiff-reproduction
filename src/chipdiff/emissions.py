'''
Emission-weight calculation for the ChIPDiff HMM
'''
import numpy as np
from scipy.integrate import quad
from scipy.special import betaln
from scipy.stats import beta as beta_distribution

from .intensity import ALPHA,TAU

def calculate_state_masses(esc_alpha, esc_beta, npc_alpha, npc_beta, tau=TAU):
    '''
    Calculate the probability masses of the joint ESC/NPC intensity distribution 
    within each ChIPDiff state. 
    States:
        non-differential: 1/tau <= p_ESC / p_NPC <= tau
        ESC-enriched:     p_ESC / p_NPC > tau
        NPC-enriched:     p_ESC / p_NPC < 1/tau

    A transformed integration variable is used for numerical stability
    because the Beta distributions are highly concentrated near zero.
    '''

    def get_npc_intensity(z):
        return -np.expm1(-z/npc_beta)
    
    def npc_transformed_density(z):
        if z == 0:
            if npc_alpha > 1:
                return 0.0
            
            return np.exp(-betaln(npc_alpha,npc_beta)-np.log(npc_beta))

        p_npc = get_npc_intensity(z)
        log_density = ((npc_alpha -1)*np.log(p_npc)-z-betaln(npc_alpha,npc_beta)-np.log(npc_beta))
        return np.exp(log_density)
    
    def non_differential_integrand(z):
        p_npc = get_npc_intensity(z)
        lower_threshold = p_npc/tau
        upper_threshold = tau*p_npc
        probability_esc_between_thresholds = (beta_distribution.cdf(upper_threshold,esc_alpha,esc_beta)
                                             - beta_distribution.cdf(lower_threshold,esc_alpha,esc_beta))
        return npc_transformed_density(z)*probability_esc_between_thresholds
    
    def esc_enriched_integrand(z):
        p_npc = get_npc_intensity(z)
        esc_threshold = tau*p_npc
        probability_esc_above_threshold = beta_distribution.sf(esc_threshold,esc_alpha,esc_beta)
        return npc_transformed_density(z)*probability_esc_above_threshold

    def npc_enriched_integrand(z):
        p_npc = get_npc_intensity(z)
        npc_threshold = p_npc/tau
        probability_esc_below_threshold = beta_distribution.cdf(npc_threshold,esc_alpha,esc_beta)
        return npc_transformed_density(z)*probability_esc_below_threshold

    non_differential_mass = quad(non_differential_integrand,0,np.inf,epsabs=1e-12,epsrel=1e-10,limit=200)[0]
    esc_enriched_mass = quad(esc_enriched_integrand,0,np.inf,epsabs=1e-12,epsrel=1e-10,limit=200)[0]
    npc_enriched_mass = quad(npc_enriched_integrand,0,np.inf,epsabs=1e-12,epsrel=1e-10,limit=200)[0]

    return(non_differential_mass,esc_enriched_mass,npc_enriched_mass)


def calculate_prior_state_masses(total_bins):
    '''
    Calculate the prior probability mass associated with each
    ChIPDiff state using the Beta(1,m) prior.
    '''
    prior_state_masses = np.array(calculate_state_masses(ALPHA,total_bins,ALPHA,total_bins,TAU))
    return prior_state_masses

def build_emission_lookup(putative_sites, total_bins):
    """
    Build a lookup table of relative HMM emission weights for each
    unique ESC/NPC fragment-count combination.

    Emission weights are calculated as posterior state mass divided
    by prior state mass.

    A state-independent observation factor is omitted because it
    cancels during HMM inference.
    """

    # Keep one row for each unique ESC/NPC count combination.
    emission_lookup = (
        putative_sites[
            [
                "esc_counts",
                "npc_counts",
                "esc_alpha_post",
                "esc_beta_post",
                "npc_alpha_post",
                "npc_beta_post"
            ]
        ]
        .drop_duplicates(
            subset=[
                "esc_counts",
                "npc_counts"
            ]
        )
        .reset_index(drop=True)
    )

    # Calculate the probability mass of each state under the prior.
    prior_state_masses = calculate_prior_state_masses(total_bins)

    emission_weights = []

    # Calculate posterior state masses for every unique count pair.
    for row in emission_lookup.itertuples(index=False):

        posterior_state_masses = np.array(calculate_state_masses(
        row.esc_alpha_post,
        row.esc_beta_post,
        row.npc_alpha_post,
        row.npc_beta_post)
        )

        # Convert posterior state masses into relative emission weights.
        relative_emissions = posterior_state_masses/prior_state_masses

        emission_weights.append(relative_emissions)

    # Convert the list of 3-element vectors into an N x 3 array.
    emission_weights = np.array(emission_weights)

    # Add the three HMM emission weights to the lookup table.
    emission_lookup[["non_diff_emission","esc_emission","npc_emission"]] = emission_weights

    return emission_lookup[["esc_counts","npc_counts","non_diff_emission","esc_emission","npc_emission"]]

def add_emission_weights(putative_sites,emission_lookup):
    '''
    Add precomputed emission weights to every putative genomic bin
    '''
    putative_sites = putative_sites.merge(emission_lookup,on=['esc_counts','npc_counts'],how='left')
    return putative_sites

def build_emission_table(putative_sites,total_bins):
    '''
    Build the emission lookup table and attach emission weights
    to all putative modification sites.
    '''
    emission_lookup = build_emission_lookup(putative_sites,total_bins)
    putative_sites = add_emission_weights(putative_sites,emission_lookup)

    return putative_sites