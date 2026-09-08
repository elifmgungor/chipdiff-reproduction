'''
End-to-end ChIPDiff reproduction pipeline.
This module connects all analysis stages:
1- raw aligned reads
2- preprocessing and genomic binning
3- putative-site filtering
4- Bayesian intensity estimation
5- HMM emission calculation
6- HMM training and inference
7- final DHMS regions
'''

from .preprocessing import preprocess_library, count_fragment_per_bin, merge_library_counts
from .filtering import load_chromosome_sizes, calculate_total_genomic_bins, build_putative_regions
from .intensity import build_intensity_table
from .emissions import build_emission_table
from .hmm import build_hmm_table
from .results import build_dhms_regions

def run_chipdiff_pipeline(esc_file,npc_file,chromosome_sizes_file, n_training_regions = 10000, random_state = 42, tolerance = 1e-6, max_iterations = 100, rho = 0.95, verbose = True):
    '''
    Run the complete ChIPdiff reproduction pipeline.
    Parameters:
        - esc_file: Path to the aligned ESC H3K27me3 reads.
        - npc_file : Path to the aligned NPC H3K27me3 reads.
        - chromosome_sizes_file: Path to the chromosome-size reference file.
        - n_training_regions: Number of randomly selected putative regions
        used for Baum-Welch HMM training.
        - random_state: Random seed used when selecting training regions.
        - tolerance: Baum-Welch convergence threshold.
        - max_iterations: Maximum number of Baum-Welch iterations.
        - rho: Posterior-probability threshold used for final DHMS calls.
        - verbose: Whether to print progress information during HMM training

    Returns:
        -results: Dictionnary containing intermediate and final analysis tables,
        HMM parameters, and useful summary values.
    '''
    # 1. Preprocess the two ChIP-seq libraries

    if verbose:
        print("Preprocessing ESC reads...")
    esc_reads = preprocess_library(esc_file)
    if verbose:
        print("Preprocessing NPC reads...")
    npc_reads = preprocess_library(npc_file)

    # Number of retained fragments afer preprocessing
    # These values are used for sequencing-depth normalization and Bayesian intensity estimation.

    n_esc = len(esc_reads)
    n_npc = len(npc_reads)

    if verbose:
        print(f'ESC fragments: {n_esc:,}')
        print(f'NPC fragments: {n_npc:,}')

    # 2. Count fragment centers in 1kb genomic bins
    esc_counts = count_fragment_per_bin(esc_reads,"esc_counts")
    npc_counts = count_fragment_per_bin(npc_reads,"npc_counts")
    bin_counts = merge_library_counts(esc_counts, npc_counts)
    if verbose:
        print(f'Observed genomic bins:{len(bin_counts):,}')

    #3. Load chromosome sizes and calculate total genomic bins.
    chromosome_sizes = load_chromosome_sizes(chromosome_sizes_file)

    # Keep the same 21 nuclear chromosomes used in the notebook analysis
    canonical_chromosomes = [f'chr{i}' for i in range(1,20)] + ['chrX','chrY']

    chromosome_sizes = chromosome_sizes[
        chromosome_sizes['chromosome'].isin(canonical_chromosomes)
    ].copy().reset_index(drop=True)

    total_bins = calculate_total_genomic_bins(chromosome_sizes)
    if verbose:
        print(f"Total genomic bins {total_bins:,}")

    #4. Identify putative histone-modification regions

    putative_sites = build_putative_regions(bin_counts,chromosome_sizes,n_esc,n_npc)
    if verbose:
        print(f"Putative bins retained {len(putative_sites):,}")
        print(f"Putative regions {putative_sites['region_id'].nunique():,}")

    #5. Estimate ESC and NPC modification intensities
    intensity_table = build_intensity_table(putative_sites,total_bins,n_esc,n_npc)

    #6 Calculate HMM emission weights
    if verbose:
        print("Calculating HMM emission weights...")

    emission_table = build_emission_table(intensity_table,total_bins)

    #7. Train the HMM infer posterior state probabilities
    if verbose:
        print("Training HMM transition matrix...")

    hmm_table, trained_transition_matrix, transition_history, training_region_ids = build_hmm_table(emission_table, n_training_regions = n_training_regions, random_state=random_state, tolerance=tolerance, max_iterations = max_iterations, rho=rho, verbose=verbose)

    #8. Build final differential histone-modification regions

    dhms_bins, dhms_regions = build_dhms_regions(hmm_table)

    if verbose:
        print()
        print("Final ChIPDiff results")
        print("----------------------")
        print(f"Differential bins {len(dhms_bins):,}")
        print(f"DHMS regions: {len(dhms_regions):,}")
        if not dhms_regions.empty:
            state_counts = dhms_regions['dhms_state'].value_counts()
            print(f"ESC-enriched regions: {state_counts.get('ESC_enriched',0):,}")
            print(f"NPC-enriched regions: {state_counts.get('NPC_enriched',0):,}")

    #9. Return intermediate and final outputs

    results = ({"bin_counts": bin_counts,
                'putative_sites':putative_sites,
                "intensity_table":intensity_table,
                "emission_table":emission_table,
                "hmm_table":hmm_table,
                "dhms_bins":dhms_bins,
                "dhms_regions":dhms_regions,
                "trained_transition_matrix":trained_transition_matrix,
                "transition_history": transition_history,
                "training_region_ids":training_region_ids,
                "n_esc":n_esc,
                "n_npc":n_npc,
                "total_bins":total_bins} )
    return results