# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    #create a list of species alignment objects, and use the list.sort() method to sort them by alignment score
    hs_gg = NeedlemanWunsch('../substitution_matrices/BLOSUM62.mat', -10, -1).align(hs_seq, gg_seq)
    hs_mm = NeedlemanWunsch('../substitution_matrices/BLOSUM62.mat', -10, -1).align(hs_seq, mm_seq)
    hs_br = NeedlemanWunsch('../substitution_matrices/BLOSUM62.mat', -10, -1).align(hs_seq, br_seq)
    hs_tt = NeedlemanWunsch('../substitution_matrices/BLOSUM62.mat', -10, -1).align(hs_seq, tt_seq)
    species_alignments = [hs_gg, hs_mm, hs_br, hs_tt] #list of objects
    species_alignment.sort(key = lambda x: x.alignment_score, reverse=True) #sort in ascending order
    d = {:"Gallus gallus", "Mus musculus", :"Balaeniceps rex", "Tursiops truncatus"}
    print("Species in order of similarity to humans (most->least):")
    print(species_alignment) #prints names of species taken from dictionary

    # print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    print("Alignment score between humans and Gallus gallus: " + str(hs_gg.alignment_score))
    print("Alignment score between humans and Mus musculus: " + str(hs_mm.alignment_score))
    print("Alignment score between humans and Balaeniceps rex: " + str(hs_br.alignment_score))
    print("Alignment score between humans and Trusiops truncatus: " + str(hs_tt.alignment_score))

if __name__ == "__main__":
    main()
