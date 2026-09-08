'''
Hidden Markov Model implementation for the ChIPDiff Reproduction
The HMM contains three hidden states:
    0. non-differential
    1. ESC-enriched
    2. NPC-enriched

Emission weights are calculated beforehand in emissions.py

The initial hidden state S0 is fixed to the non-differential state.
Therefore, the first genomic bin is reached through a transition from 
state 0 rather than being forced to state 0 itself.
'''
import numpy as np

STATES =['non_differential','ESC_enriched','NPC_enriched']
EMISSION_COLUMNS = ['non_diff_emission','esc_emission','npc_emission']

RHO = 0.95

def initialize_transition_matrix():
    '''
    Initialize the 3x3 HMM transition matrix uniformly.
    ChIPDiff initializes every transition probability to 1/3
    '''
    return np.full((3,3),1/3,dtype=float)

def build_training_sequences(putative_sites,n_regions=10000,random_state=42):
    '''
    Randomly select putative modification regions and convert them into HMM
    emission sequences.
    
    Parameters
    - putative_sites: DataFrame containing region IDs and HMM emission weights
    - n_regions: Number of putative regions used for Baum-Welch training.
    The ChIPDiff paper uses 10,000 regions.
    - random_state: Random seed used to make the sampling reproducible

    Returns
    -training_sequences: List of arrays with shape: number_of_bins_in_region x3 
    where the three columns correspond to three HMM states
    -training_region_ids: IDs of the regions selected for training.
    '''

    #Obtain one entry per putative modification region
    all_region_ids = putative_sites['region_id'].drop_duplicates()

    # Avoid requesting more regions than are available
    n_regions = min(n_regions, len(all_region_ids))

    #Randomly sample regions for Baum-Welch training
    training_region_ids = all_region_ids.sample(n=n_regions,random_state=random_state).to_numpy()

    #keep only bins belonging to selected regions
    training_bins = putative_sites[putative_sites['region_id'].isin(training_region_ids)].sort_values(['region_id','chromosome','bin'])

    training_sequences =[]

    #Each putative region becomes one independent HMM sequence
    for _, region in training_bins.groupby("region_id",sort=False):
        emission_sequence = region[EMISSION_COLUMNS].to_numpy(dtype=float)
        training_sequences.append(emission_sequence)

    return(training_sequences, training_region_ids)

def forward_pass(emission_sequence,transition_matrix):
    '''
    Perform the scaled HMM forward algorithm for one genomic region.
    The state before the first genomic bin (S0) is fixed to non differential state. 
    Therefore, probabilities for the first genomic bin are obtained from row 0 of the transition matrix
    '''
    n_bins = emission_sequence.shape[0]
    n_states = transition_matrix.shape[0]

    forward_probs = np.zeros((n_bins,n_states))
    scaling_factors = np.zeros(n_bins)

    #S0 is fixed to state 0. The first genomic bin is reached through the 
    #transitions from S0

    first_state_probs = transition_matrix[0]
    forward_probs[0] = first_state_probs * emission_sequence[0]
    scaling_factors[0] = forward_probs[0].sum()

    if scaling_factors[0] <= 0:
        raise ValueError("Invalid forward scaling factor at the first genomic bin")
    
    forward_probs[0] /= scaling_factors[0]

    # Propagate probabilities through the remaining genomic bins
    for t in range(1,n_bins):
        predicted_state_probs = forward_probs[t-1] @ transition_matrix
        forward_probs[t] = predicted_state_probs * emission_sequence[t]
        scaling_factors[t] = forward_probs[t].sum()

        if scaling_factors[t] <= 0:
            raise ValueError(f'Invalid forward scaling factor at genomic bin index {t}')

        forward_probs[t] /= scaling_factors[t]

    return(forward_probs, scaling_factors)

def backward_pass(emission_sequence, transition_matrix, scaling_factors):
    '''
    Perform the scaled HMM backward algorithm for one genomic region
    '''
    n_bins = emission_sequence.shape[0]
    n_states = transition_matrix.shape[0]

    backward_probs = np.zeros((n_bins,n_states))
    
    # No observations remain after the final genomic bin
    backward_probs[-1] = 1.0

    # Move backward through the sequence
    for t in range (n_bins -2, -1,-1):

        next_bin_information = emission_sequence[t+1] * backward_probs[t+1]   
        backward_probs[t] = transition_matrix @ next_bin_information

        # Use the same scaling factors as the forward algorithm
        backward_probs[t] /= scaling_factors[t+1]
    
    return backward_probs

def calculate_posterior_probabilities(forward_probs, backward_probs):
    '''
    Calculate posterior HMM state probabilities
    '''
    posterior_probs = forward_probs * backward_probs

    # Normalize each genomic bin so that its three state probabilities sum to one.

    posterior_probs /= posterior_probs.sum(axis=1,keepdims=True)
    return posterior_probs

def calculate_xi(emission_sequence, transition_matrix, forward_probs, backward_probs):
    '''
    Calculate expected transitions between neighboring genomic bins.
    xi[t,i,j] represents the posterior probability that:
    bin t is in state i and bin t+1 is in state j
    '''
    n_bins = emission_sequence.shape[0]
    n_states = transition_matrix.shape[0]

    # A region containing n bins contains n-1 transitions between neighboring genomic bins
    xi = np.zeros((n_bins -1, n_states, n_states))

    for t in range(n_bins -1):
        transition_scores = (forward_probs[t][:,np.newaxis] 
                            * transition_matrix 
                            * emission_sequence[t+1][np.newaxis,:]
                            * backward_probs[t+1][np.newaxis,:]
                            )
        total_score = transition_scores.sum()
        if total_score <=0:
            raise ValueError(f'Could not normalize xi at transition index {t}')
        
        xi[t] = transition_scores / total_score
    
    return xi

def update_transition_matrix(training_sequences,transition_matrix):
    '''
    Perform one Baum-Welch update of the transition matrix
    Expected transition counts are accumulated across all training regions
    before the matrix is row-normalized
    '''
    n_states = transition_matrix.shape[0]

    transition_counts = np.zeros((n_states,n_states))

    for emission_sequence in training_sequences :

        forward_probs, scaling_factors = forward_pass(emission_sequence,transition_matrix)
        backward_probs = backward_pass(emission_sequence,transition_matrix,scaling_factors)
        posterior_probs = calculate_posterior_probabilities(forward_probs,backward_probs)

        #Every putative region starts from the fixed S0 state, which is non differential state (state 0)
        # Therefore gamma for the first bin contributes expected transitions from S0 to each possible first bin state

        transition_counts[0] += posterior_probs[0]
        
        #Regions containing more than one bin also contribute to transitions between neighboring bins
        if emission_sequence.shape[0] >1:
            xi = calculate_xi(emission_sequence, transition_matrix, forward_probs, backward_probs)
            transition_counts += xi.sum(axis=0)
    
    row_sums = transition_counts.sum(axis=1, keepdims=True)

    if np.any(row_sums == 0):
        raise ValueError('At least one HMM state received no expected outgoing transitions')
    
    updated_transition_matrix = transition_counts / row_sums

    return updated_transition_matrix

def train_transition_matrix(training_sequences,initial_transition_matrix=None,tolerance=1e-6,max_iterations=100,verbose=True):
    '''
    Train the HMM transition matrix using Baum-Welch iterations.
    Training stops when the largest absolute change in any transition probability 
    falls below the convergence tolerance.

    The tolerance and maximum number of iterations are implementation choices for this reproduction
    because the original paper does not report explicit values for them
    '''
    if initial_transition_matrix is None:
        transition_matrix = initialize_transition_matrix()
    else:
        transition_matrix = initial_transition_matrix.copy().astype(float)
    
    transition_history = []
    
    for iteration in range(max_iterations):
        
        updated_transition_matrix = update_transition_matrix(training_sequences,transition_matrix)
        
        #Largest change among the nine transition probabilities
        max_change = np.max(np.abs(updated_transition_matrix - transition_matrix))
        transition_history.append(max_change)

        transition_matrix = updated_transition_matrix

        if verbose:
            print(f'Iteration {iteration+1} max change {max_change:.8f}')

        if max_change < tolerance:
            if verbose:
                print(f'Baum-Welch converged after {iteration+1} iterations')
            break
    
    return(transition_matrix,transition_history)

def infer_hmm_states(putative_sites,transition_matrix,rho=RHO):
    '''
    Perform final HMM inference for all putative modification regions.

    Forward-backward inference is performed independently within each putative region.

    A genomic bin is called ESC-enriched or NPC-enriched only when
    a corresponding posterior state probability exceeds rho.
    '''
    
    result = putative_sites.sort_values(['region_id','chromosome','bin']).reset_index(drop=True).copy()   
    n_bins = len(result)
    posterior_matrix = np.zeros((n_bins,3))
    
    #Run forward-backward separately for every putative region
    for _, region in result.groupby("region_id",sort=False):
        
        row_indices = region.index.to_numpy()
        emission_sequence = region[EMISSION_COLUMNS].to_numpy(dtype=float)

        forward_probs, scaling_factors = forward_pass(emission_sequence,transition_matrix)
        backward_probs = backward_pass(emission_sequence, transition_matrix, scaling_factors)
        posterior_probs = calculate_posterior_probabilities(forward_probs,backward_probs)

        posterior_matrix[row_indices] = posterior_probs
    
    # Store posterior probabilities as explicit DataFrame columns.
    
    result['non_diff_posterior'] = posterior_matrix[:,0]
    result['esc_posterior'] = posterior_matrix[:,1]
    result['npc_posterior'] = posterior_matrix[:,2]

    # Maximum posterior state, without applying the DHMS threshold
    state_indices = np.argmax(posterior_matrix, axis=1)
    result['hmm_state'] = np.array(STATES)[state_indices]

    #Final ChIPDiff differential calls. 
    # A bin must have posterior probability > rho for one of the differential states
    # Otherwise it remains non-differential

    result['dhms_state'] = "non_differential"
    result.loc[result['esc_posterior']>rho,'dhms_state'] = 'ESC_enriched'
    result.loc[result['npc_posterior']>rho,'dhms_state'] = 'NPC_enriched'

    return result

def build_hmm_table(putative_sites, n_training_regions=10000, random_state=42, tolerance=1e-6, max_iterations=100, rho=RHO, verbose=True):
    '''
    Run HMM training and final inference.
    This is the high-level function for the HMM stage of the ChIPDiff reproduction
    '''

    training_sequences, training_region_ids = build_training_sequences(putative_sites,n_regions=n_training_regions,random_state=random_state)
    initial_transition_matrix = initialize_transition_matrix()
    trained_transition_matrix, transition_history = train_transition_matrix(training_sequences,initial_transition_matrix=initial_transition_matrix,tolerance=tolerance,max_iterations=max_iterations,verbose=verbose)
    hmm_table = infer_hmm_states(putative_sites,trained_transition_matrix,rho=rho)

    return(hmm_table, trained_transition_matrix, transition_history, training_region_ids)