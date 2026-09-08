'''
Command-line entry point for the ChIPDiff reproduction pipeline
'''
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(project_root/"src"))
from chipdiff.pipeline import run_chipdiff_pipeline

def main():
    esc_file = project_root/'data/raw/GSM307619_ES.H3K27me3.aligned.txt.gz'
    npc_file = project_root/'data/raw/GSM307614_NP.H3K27me3.aligned.txt.gz'
    chromosome_sizes_file = project_root/'metadata/mm8.chrom.sizes'

    results =  run_chipdiff_pipeline(esc_file = esc_file, npc_file=npc_file,
               chromosome_sizes_file = chromosome_sizes_file,
               n_training_regions = 10000, random_state = 42,
               tolerance = 1e-6, max_iterations=100,
               rho = 0.95, verbose=True)
    
    results_directory = project_root/'results/tables'
    results_directory.mkdir(parents=True,exist_ok=True)

    #Save final differetial genomic bins
    results['dhms_bins'].to_csv(results_directory/"dhms_bins.csv",index=False)

    #Save final merged DHMS regions
    results['dhms_regions'].to_csv(results_directory/"dhms_regions.csv",index=False)

    #Save the complete HMM bin-level output
    results['hmm_table'].to_csv(results_directory/'hmm_results.csv',index=False)

    print()
    print("Results saved to",results_directory)

if __name__ == "__main__":
    main()