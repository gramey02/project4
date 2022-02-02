# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    #test for proper matrix filling by assessing if the filled matrices equal what is expected
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    #create a NW object
    NW = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10, -1)
    #get alignment between seq1 and seq2 using the align function
    NW.align(seq1, seq2)
    
    #check that the alignment matrix equals what it should
    expected = np.array([[  0., -inf, -inf, -inf],
                         [-inf,   5., -11., -13.],
                         [-inf, -12.,   4.,  -8.],
                         [-inf, -12.,  -1.,   5.],
                         [-inf, -14.,  -6.,   4.]])
    observed = NW._align_matrix
    comparison = observed==expected
    observed_equals_expected = comparison.all()
    assert observed_equals_expected == True
    
    #check that the gap_B matrix equals what it should
    expected = np.array([[-10., -inf, -inf, -inf],
                         [-11., -22., -23., -24.],
                         [-12.,  -6., -17., -18.],
                         [-13.,  -7.,  -7., -18.],
                         [-14.,  -8.,  -8.,  -6.]])
    #check that the gap_A matrix equals what it should
    
    #test a few other alignment outputs, just to see that everything is working is it should
    assert NW.alignment_score == 4 #alignment score should be 4
    assert NW.seqA_align == "MYQR" #make sure seqA alignment is correct
    assert NW.seqB_align == "M-QR" #make sure seqB alignment is correct
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    pass




