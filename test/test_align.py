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
    observed = NW._gapB_matrix
    comparison = observed==expected
    observed_equals_expected = comparison.all()
    assert observed_equals_expected == True
    
    #check that the gap_A matrix equals what it should
    expected = np.array([[-10., -11., -12., -13.],
                         [-inf, -22.,  -6.,  -7.],
                         [-inf, -23., -17.,  -7.],
                         [-inf, -24., -18., -12.],
                         [-inf, -25., -19., -17.]])
    observed = NW._gapA_matrix
    comparison = observed==expected
    observed_equals_expected = comparison.all()
    assert observed_equals_expected == True
    
    #test a few other alignment outputs, just to check that everything is working is it should be
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
    #test that the backtrace matrices are filled properly
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    NW2 = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10, -1)
    NW2.align(seq3, seq4)
    
    #check that each of the backtrace matrices equals what it should
    #_back check
    expected = np.array([[('align_matrix', 0, 0), ('align_matrix', 0, 1), ('align_matrix', 0, 2), -inf],
                         [('align_matrix', 1, 0), ('align_matrix', 0, 0), ('gapA_matrix', 0, 1), ('gapA_matrix', 0, 2)],
                         [('align_matrix', 2, 0), ('gapB_matrix', 1, 0), ('align_matrix', 1, 1), ('gapA_matrix', 1, 2)],
                         [('align_matrix', 3, 0), ('gapB_matrix', 2, 0), ('gapB_matrix', 2, 1), ('align_matrix', 2, 2)],
                         [-inf, ('gapB_matrix', 3, 0), ('gapB_matrix', 3, 1), ('align_matrix', 3, 2)]], dtype=object)
    observed = NW._back
    comparison = observed==expected
    observed_equals_expected = comparison.all()
    assert observed_equals_expected == True
    
    #_back_B check
    expected = np.array([[('gapB_matrix', 0, 0), ('gapB_matrix', 0, 1), ('gapB_matrix', 0, 2), ('gapB_matrix', 0, 3)],
                         [('gapB_matrix', 1, 0), ('gapA_matrix', 0, 1), ('gapA_matrix', 0, 2), ('gapA_matrix', 0, 3)],
                         [('gapB_matrix', 2, 0), ('align_matrix', 1, 1), ('gapA_matrix', 1, 2), ('gapA_matrix', 1, 3)],
                         [('gapB_matrix', 3, 0), ('gapB_matrix', 2, 1), ('align_matrix', 2, 2), ('gapA_matrix', 2, 3)],
                         [('gapB_matrix', 4, 0), ('gapB_matrix', 3, 1), ('gapB_matrix', 3, 2), ('align_matrix', 3, 3)]],
                        dtype=object)
    observed = NW._back_B
    comparison = observed==expected
    observed_equals_expected = comparison.all()
    assert observed_equals_expected == True
    
    #_back_A check
    expected = np.array([[('gapA_matrix', 0, 0), ('gapA_matrix', 0, 1), ('gapA_matrix', 0, 2), ('gapA_matrix', 0, 3)],
                         [('gapA_matrix', 1, 0), ('gapB_matrix', 1, 0), ('align_matrix', 1, 1), ('gapA_matrix', 1, 2)],
                         [('gapA_matrix', 2, 0), ('gapB_matrix', 2, 0), ('gapB_matrix', 2, 1), ('align_matrix', 2, 2)],
                         [('gapA_matrix', 3, 0), ('gapB_matrix', 3, 0), ('gapB_matrix', 3, 1), ('align_matrix', 3, 2)],
                         [('gapA_matrix', 4, 0), ('gapB_matrix', 4, 0), ('gapB_matrix', 4, 1), ('align_matrix', 4, 2)]],
                        dtype=object)
    observed = NW._back_A
    comparison = observed==expected
    observed_equals_expected = comparison.all()
    assert observed_equals_expected == True
    
    




