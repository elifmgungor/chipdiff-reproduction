from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]

RESULTS_DIRECTORY = PROJECT_ROOT /'results/tables'
FIGURES_DIRECTORY = PROJECT_ROOT /'results/figures'

DHMS_REGIONS_FILE = RESULTS_DIRECTORY / 'dhms_regions.csv'
DHMS_BINS_FILE = RESULTS_DIRECTORY / 'dhms_bins.csv'
ANNOTATED_REGIONS_FILE = RESULTS_DIRECTORY / 'dhms_regions_annotated.csv'


def load_results():
    dhms_bins = pd.read_csv(DHMS_BINS_FILE)
    dhms_regions = pd.read_csv(DHMS_REGIONS_FILE)
    annotated_regions = pd.read_csv(ANNOTATED_REGIONS_FILE)
    return dhms_bins, dhms_regions, annotated_regions


def plot_region_length_distribution(dhms_regions):
    """
    Plot the length distribution of reproduced H3K27me3 DHMS regions
    using the categories reported in the original ChIPDiff study.
    """
    region_length_kb = ( dhms_regions["region_end"] - dhms_regions["region_start"] ) / 1000

    counts = [(region_length_kb == 1).sum(), (region_length_kb == 2).sum(),
    ((region_length_kb > 2) & (region_length_kb <= 5)).sum(),
    ( (region_length_kb > 5) & (region_length_kb <= 10)).sum(),
    (region_length_kb > 10).sum() ]

    labels = [ "1 kb","2 kb","3-5 kb","6-10 kb",">10 kb"]

    fig, ax = plt.subplots(figsize=(8,5))

    ax.bar(labels, counts)
    ax.set_ylabel("DHMS regions")
    ax.set_title("H3K27me3 DHMS region-length distribution")
    fig.tight_layout()
    fig.savefig( FIGURES_DIRECTORY / "region_length_distribution.png",dpi=300)
    plt.close(fig)


def plot_intensity_scatter(dhms_bins):
    """
    Plot posterior mean H3K27me3 intensities in ESC and NPC
    for differential bins identified by the HMM.
    """
    fig, ax = plt.subplots(figsize=(6.5,6.5))

    for state, state_bins in dhms_bins.groupby("dhms_state"):

        ax.scatter(np.log10(state_bins["esc_intensity"]),
            np.log10(state_bins["npc_intensity"]),
            s=7,alpha=0.35,label=state
        )

    minimum = min(np.log10(dhms_bins["esc_intensity"]).min(), 
                  np.log10(dhms_bins["npc_intensity"]).min())

    maximum = max(
        np.log10(dhms_bins["esc_intensity"]).max(),
        np.log10(dhms_bins["npc_intensity"]).max())

    ax.plot([minimum,maximum],[minimum,maximum], linestyle="--", linewidth=1)

    ax.set_xlim(minimum,maximum)
    ax.set_ylim(minimum,maximum)

    ax.set_xlabel( "log10 ESC posterior mean intensity")
    ax.set_ylabel("log10 NPC posterior mean intensity")

    ax.set_title("Differential H3K27me3 bins")
    ax.legend()

    fig.tight_layout()
    fig.savefig(FIGURES_DIRECTORY/ "intensity_scatter.png",dpi=300)
    plt.close(fig)


def plot_genomic_annotation(annotated_regions):
    """
    Plot genomic context of ESC- and NPC-enriched DHMS regions.
    """
    annotation_order = ["promoter","gene_body","intergenic"]
    state_order = ["ESC_enriched","NPC_enriched"]
    annotation_counts = pd.crosstab( annotated_regions["dhms_state"], annotated_regions["annotation"])
    annotation_counts = annotation_counts.reindex(index=state_order, 
    columns=annotation_order, fill_value=0)

    annotation_percentages = (annotation_counts.div(annotation_counts.sum(axis=1),axis=0)* 100)
    labels = ["ESC-enriched","NPC-enriched"]

    x = np.arange(len(labels))
    bottom = np.zeros(len(labels))
    fig, ax = plt.subplots(figsize=(7,5))

    for annotation in annotation_order:

        values = annotation_percentages[annotation].to_numpy()
        ax.bar(x,values,bottom=bottom,label=annotation)
        bottom += values

    ax.set_ylabel("DHMS regions (%)")
    ax.set_title("Genomic context of differential H3K27me3 regions")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0,100)
    ax.legend(title="Annotation")
    fig.tight_layout()
    fig.savefig( FIGURES_DIRECTORY / "genomic_annotation.png",dpi=300)
    plt.close(fig)


def main():
    FIGURES_DIRECTORY.mkdir(parents=True, exist_ok=True)
    dhms_bins, dhms_regions, annotated_regions = load_results()
    plot_region_length_distribution(dhms_regions)
    plot_intensity_scatter(dhms_bins)
    plot_genomic_annotation(annotated_regions)
    print(f"Figures saved to {FIGURES_DIRECTORY}")

if __name__ == "__main__":
    main()