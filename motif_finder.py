import statistics
import math

class Motif:
    """
    The counts parameter is a dictionary, that maps each nucleotide to a list of integers. Each integer indicates the number of times that the nucleotide appeared in the column.
    
    background_frequencies is a dictionary that maps the nucleotide to their frequence in the dataset.

    num_sequences is the number of sequences that went into creating the counts


    pseudocount_value is the value used for the pseudocounts (recommend 0.01)
    """
    def __init__(self, sequences, background_frequencies, pseudocount_value = 0.01):
        self.background_frequencies = background_frequencies
        amino_acids = list(background_frequencies.keys())
        pwm = dict()
        num_sequences = len(sequences)
        width = len(sequences[0])
        denom = pseudocount_value  + num_sequences
        for position in range(0, width):
            probs = {amino_acid: 0.0 for amino_acid in amino_acids}
            for seq in sequences:
                probs[seq[position]] += 1.0

            for aa in amino_acids:
                residue_prob_background = background_frequences[aa]
                pwm[aa].append(math.log((probs[aa] + residue_prob_background*pseudocount_value)/(denom*residue_prob_background)))
        self.pwm = pwm
        self.pwm_width = width
        self.pseudocount_value = pseudocount_value

        
    def get_pwm(self):
        return self.pwm
    def recalculate(self, sequences):
        self = Motif(sequences, self.background_frequencies, self.pseudocount_value)


    @staticmethod
    def compute_background_frequencies(genome):
        #genome_only_nuc = list(filter(lambda x: x in ['A', 'G', 'C', 'T'], genome))
        genome_size = len(genome)
        return dict(map(lambda x: (x[0], 1.0*x[1]/genome_size), Counter(genome).items()))

    def get_pwm_width(self):
        return self.pwm_width

    def likelihood(self, sequence):
        likelihood = 0
        i = 0
        for x in sequence:
           likelihood += self.pwm[x][i]
           i += 1
        return math.exp(likelihood)


    

def cluster_motifs(sequences, motifs):
    """
    sequences is a list of sequences
    motifs is a list of motifs (each of which is an instance of the Motif class)
    This function goes through the sequences, and for each sequence, finds the motif that gives the highest likelihood. 
    It returns the clusters.u and likelihood
    TODO: Worry about possible multiple hypothesis testing (this should help us constrain the number of motifs)
    """
    clusters = {motif: [] for motif in motifs}
    likelihoods = []
    for seq in sequences:
        highest_likelihood = motifs[0].likelihood(seq)
        best_motif = motifs[0]
        for motif in motifs[1::]:
            likelihood = motif.likelihood(seq)
            if likelihood > highest_likelihood:
                highest_likelihood = likelihood
                best_motif = motif
        clusters[best_motif].append(seq)
        likelihoods.append(highest_likelihood)
    
    return {'clusters': clusters, 'likelihood': statistics.mean(likelihoods)}

def compute_motifs(clusters):
    #for each motif, re-compute the parameters
    for motif, sequences in clusters:
        motif.recalculate(sequences)

#compute the background frequencies
def background_frequencies(proteome_file):
    
