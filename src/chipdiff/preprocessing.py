'''
Preprocessing utilities for the ChIPDiff reproduction.
This module implements the preprocessing steps applied to aligned single-end Chip-seq reads before putative histone-modification-site detection
'''

from pathlib import Path
import numpy as np
import pandas as pd

ALIGNED_COLUMNS = ['chromosome','start','end','orientation','read_id','mismatches','sequence']
DEFAULT_FRAGMENT_SHIFT = 100
DEFAULT_BIN_SIZE = 1000
BIN_COLUMNS = ["chromosome","bin","bin_start","bin_end"]

def load_aligned_reads(file_path):
    '''
    Load an aligned. ChIP seq-read file into a pandas DataFrame
    The aligned files containes seven tab-seperated columns:
    chromosome, start, end, orientation, read_id, mismatches and sequence
    '''
    reads = pd.read_csv(file_path, sep='\t', names= ALIGNED_COLUMNS)
    return reads

def define_tag_position(df):
    '''
    Define genomic tag position according to strand orientation.
    Forward tags use the alignment start position
    Backward tags use end-1 under the half-open coordinate convention.
    '''
    df= df.copy()
    df['tag_position'] = np.where(df['orientation']== '-',df['end']-1,df['start'])
    return df

def remove_redundant_tags(df):
    '''
    Remove potentially PCR-derived redundant tags.
    Tags mapped to the same chromsome, position, and orientation are counted only once.
    '''
    df = df.drop_duplicates(subset=['chromosome','orientation','tag_position'],keep='first')
    df = df.reset_index(drop=True)
    return df

def estimate_fragment_centers(df):
    '''
    Estimate the enter of each ChIP fragment by shifting the tag position by 100 bp in the direction
    of its orientation
    '''
    df = df.copy()
    df['fragment_center'] = np.where(df['orientation'] == '+',df['tag_position'] + 100, df['tag_position'] - 100)
    return df

def exclude_chrM(df):
    '''
    Exclude mitochondrial reads before genomic binning
    This is an implementation choice because shifting fragments across
    the circular mitochrondial genome prigin is not handled here
    '''
    df = df[df['chromosome']!='chrM'].copy()
    df = df.reset_index(drop=True)
    return df

def assign_genomic_bins(df):
    '''
    Assign fragment centers to non-overlaping genomic bins
    '''
    df = df.copy()
    df['bin'] = df['fragment_center'] // 1000
    df['bin_start'] = df['bin'] * 1000
    df['bin_end'] = df['bin_start'] + 1000
    return df

def count_fragment_per_bin(df,count_column):
    '''
    count fragment centers in each genomic bin
    '''
    counts = df.groupby(['chromosome','bin','bin_start','bin_end']).size().reset_index(name = count_column)
    return counts

def merge_library_counts(esc_counts, npc_counts):
    '''
    Combine ESC and NPC fragment counts into one-bin-level table.
    Bins present in only one library are retained and assigned 
    zero fragments in the other library
    '''
    bin_counts = esc_counts.merge(npc_counts, on=['chromosome','bin','bin_start','bin_end'],how='outer')
    bin_counts[['esc_counts','npc_counts']] = bin_counts[['esc_counts','npc_counts']].fillna(0).astype('int')
    return bin_counts

def preprocess_library(file_path):
    '''
    Run read-level preprocessing for one ChIP-seq library
    '''
    reads = load_aligned_reads(file_path)
    reads = define_tag_position(reads)
    reads = remove_redundant_tags(reads)
    reads = estimate_fragment_centers(reads)
    reads = exclude_chrM(reads)
    reads = assign_genomic_bins(reads)

    return reads

def build_bin_count_table(esc_file, npc_file):
    '''
    Preprocess ESC and NPC libraries and build the combined
    bin-level fragment count table.
    '''
    esc_reads = preprocess_library(esc_file)
    npc_reads = preprocess_library(npc_file)

    esc_counts = count_fragment_per_bin(esc_reads,"esc_counts")
    npc_counts = count_fragment_per_bin(npc_reads, "npc_counts")

    bin_counts = merge_library_counts(esc_counts, npc_counts)

    return bin_counts