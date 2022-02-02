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
    hs_gg = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10, -1)
    hs_gg.align(hs_seq, gg_seq) #homo sapiens and gallus gallus alignment
    hs_mm = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10, -1)
    hs_mm.align(hs_seq, mm_seq) #homo sapiens and mus musculus alignment
    hs_br = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10, -1)
    hs_br.align(hs_seq, br_seq) #homo sapiens and balaeniceps rex alignment
    hs_tt = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10, -1)
    hs_tt.align(hs_seq, tt_seq) #homo sapiens and tursiops truncatus alignment
    
    species_alignments = [hs_gg, hs_mm, hs_br, hs_tt] #list of created objects
    species_alignments.sort(key = lambda x: x.alignment_score, reverse=True) #sort object list in ascending order
    #create dictionary of alignment scores with aligned species
    d = {3173:"Gallus gallus", 3682:"Mus musculus", 2941:"Balaeniceps rex", 3916:"Tursiops truncatus"}
    print("Species in order of similarity to humans (most->least):")
    for i in range(0,4):
        print(str(i+1) + ": " + d[species_alignments[i].alignment_score]) #prints names of species taken from dictionary

    # print all of the alignment scores between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    print("Alignment score between humans and Gallus gallus: " + str(hs_gg.alignment_score))
    print("Alignment score between humans and Mus musculus: " + str(hs_mm.alignment_score))
    print("Alignment score between humans and Balaeniceps rex: " + str(hs_br.alignment_score))
    print("Alignment score between humans and Trusiops truncatus: " + str(hs_tt.alignment_score))

if __name__ == "__main__":
    main()
