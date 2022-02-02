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
    expected = np.array([[  0., -np.inf, -np.inf, -np.inf],
                         [-np.inf,   5., -11., -13.],
                         [-np.inf, -12.,   4.,  -8.],
                         [-np.inf, -12.,  -1.,   5.],
                         [-np.inf, -14.,  -6.,   4.]], dtype=object)
    observed = NW._align_matrix
    for i in range(0,len(seq1)):
        for j in range(0,len(seq2)):
            if np.isinf(expected[i,j]) == True:
                assert np.isinf(observed[i,j]) ==True
            else:
                assert expected[i,j] == observed[i,j]
    
    #check that the gap_B matrix equals what it should
    expected = np.array([[-10., -np.inf, -np.inf, -np.inf],
                         [-11., -22., -23., -24.],
                         [-12.,  -6., -17., -18.],
                         [-13.,  -7.,  -7., -18.],
                         [-14.,  -8.,  -8.,  -6.]])
    observed = NW._gapB_matrix
    for i in range(0,len(seq1)):
        for j in range(0,len(seq2)):
            if np.isinf(expected[i,j]) == True:
                assert np.isinf(observed[i,j]) ==True
            else:
                assert expected[i,j] == observed[i,j]
    
    #check that the gap_A matrix equals what it should
    expected = np.array([[-10., -11., -12., -13.],
                         [-np.inf, -22.,  -6.,  -7.],
                         [-np.inf, -23., -17.,  -7.],
                         [-np.inf, -24., -18., -12.],
                         [-np.inf, -25., -19., -17.]])
    observed = NW._gapA_matrix
    for i in range(0,len(seq1)):
        for j in range(0,len(seq2)):
            if np.isinf(expected[i,j]) == True:
                assert np.isinf(observed[i,j])==True
            else:
                assert expected[i,j] == observed[i,j]
    
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
    expected = np.array([[('align_matrix', 0, 0), ('align_matrix', 0, 1),('align_matrix', 0, 2), ('align_matrix', 0, 3),
                          ('align_matrix', 0, 4), ('align_matrix', 0, 5),('align_matrix', 0, 6), -np.inf],
                         [('align_matrix', 1, 0), ('align_matrix', 0, 0), ('gapA_matrix', 0, 1), ('gapA_matrix', 0, 2),
                          ('gapA_matrix', 0, 3), ('gapA_matrix', 0, 4),('gapA_matrix', 0, 5), ('gapA_matrix', 0, 6)],
                         [('align_matrix', 2, 0), ('gapB_matrix', 1, 0),('align_matrix', 1, 1), ('gapA_matrix', 1, 2),
                          ('gapA_matrix', 1, 3), ('gapA_matrix', 1, 4),('gapA_matrix', 1, 5), ('gapA_matrix', 1, 6)],
                         [('align_matrix', 3, 0), ('gapB_matrix', 2, 0),('gapB_matrix', 2, 1), ('align_matrix', 2, 2),
                          ('align_matrix', 2, 3), ('align_matrix', 2, 4),('align_matrix', 2, 5), ('gapA_matrix', 2, 6)],
                         [('align_matrix', 4, 0), ('gapB_matrix', 3, 0),('gapB_matrix', 3, 1), ('gapB_matrix', 3, 2),
                          ('align_matrix', 3, 3), ('align_matrix', 3, 4),('gapA_matrix', 3, 5), ('gapA_matrix', 3, 6)],
                         [('align_matrix', 5, 0), ('gapB_matrix', 4, 0),('gapB_matrix', 4, 1), ('align_matrix', 4, 2),
                          ('gapB_matrix', 4, 3), ('align_matrix', 4, 4),('align_matrix', 4, 5), ('align_matrix', 4, 6)],
                         [('align_matrix', 6, 0), ('gapB_matrix', 5, 0),('gapB_matrix', 5, 1), ('align_matrix', 5, 2),
                          ('gapB_matrix', 5, 3), ('align_matrix', 5, 4),('align_matrix', 5, 5), ('align_matrix', 5, 6)],
                         [('align_matrix', 7, 0), ('gapB_matrix', 6, 0),('gapB_matrix', 6, 1), ('gapB_matrix', 6, 2),
                          ('align_matrix', 6, 3), ('align_matrix', 6, 4),('gapB_matrix', 6, 5), ('align_matrix', 6, 6)],
                         [('align_matrix', 8, 0), ('gapB_matrix', 7, 0),('gapB_matrix', 7, 1), ('gapB_matrix', 7, 2),
                          ('align_matrix', 7, 3), ('align_matrix', 7, 4),('gapA_matrix', 7, 5), ('gapA_matrix', 7, 6)],
                         [('align_matrix', 9, 0), ('gapB_matrix', 8, 0),('gapB_matrix', 8, 1), ('align_matrix', 8, 2),
                          ('gapB_matrix', 8, 3), ('gapB_matrix', 8, 4),('align_matrix', 8, 5), ('gapA_matrix', 8, 6)],
                         [-np.inf, ('gapB_matrix', 9, 0), ('gapB_matrix', 9, 1),('align_matrix', 9, 2), ('gapB_matrix', 9, 3),
                          ('gapB_matrix', 9, 4), ('align_matrix', 9, 5),('align_matrix', 9, 6)]], dtype=object)
    observed = NW2._back
    for i in range(0,len(seq3)):
        for j in range(0,len(seq4)):
            if (i==10 and j==0) or (i==0 and j==8):
                assert np.isninf(observed[i,j])==True
            else:
                assert expected[i,j] == observed[i,j]


    
    #_back_B check
    expected = np.array([[('gapB_matrix', 0, 0), ('gapB_matrix', 0, 1),('gapB_matrix', 0, 2), ('gapB_matrix', 0, 3),
                          ('gapB_matrix', 0, 4), ('gapB_matrix', 0, 5),('gapB_matrix', 0, 6), ('gapB_matrix', 0, 7)],
                         [('gapB_matrix', 1, 0), ('gapA_matrix', 0, 1),('gapA_matrix', 0, 2), ('gapA_matrix', 0, 3),
                          ('gapA_matrix', 0, 4), ('gapA_matrix', 0, 5),('gapA_matrix', 0, 6), ('gapA_matrix', 0, 7)],
                         [('gapB_matrix', 2, 0), ('align_matrix', 1, 1),('gapA_matrix', 1, 2), ('gapA_matrix', 1, 3),
                          ('gapA_matrix', 1, 4), ('gapA_matrix', 1, 5),('gapA_matrix', 1, 6), ('gapA_matrix', 1, 7)],
                         [('gapB_matrix', 3, 0), ('gapB_matrix', 2, 1),('align_matrix', 2, 2), ('align_matrix', 2, 3),
                          ('align_matrix', 2, 4), ('align_matrix', 2, 5),('gapA_matrix', 2, 6), ('align_matrix', 2, 7)],
                         [('gapB_matrix', 4, 0), ('gapB_matrix', 3, 1),('gapB_matrix', 3, 2), ('align_matrix', 3, 3),
                          ('align_matrix', 3, 4), ('gapA_matrix', 3, 5),('gapA_matrix', 3, 6), ('gapA_matrix', 3, 7)],
                         [('gapB_matrix', 5, 0), ('gapB_matrix', 4, 1),('gapB_matrix', 4, 2), ('gapB_matrix', 4, 3),
                          ('align_matrix', 4, 4), ('align_matrix', 4, 5),('align_matrix', 4, 6), ('align_matrix', 4, 7)],
                         [('gapB_matrix', 6, 0), ('gapB_matrix', 5, 1),('gapB_matrix', 5, 2), ('gapB_matrix', 5, 3),
                          ('gapB_matrix', 5, 4), ('align_matrix', 5, 5),('gapB_matrix', 5, 6), ('align_matrix', 5, 7)],
                         [('gapB_matrix', 7, 0), ('gapB_matrix', 6, 1),('gapB_matrix', 6, 2), ('gapB_matrix', 6, 3),
                          ('gapB_matrix', 6, 4), ('gapB_matrix', 6, 5),('align_matrix', 6, 6), ('gapB_matrix', 6, 7)],
                         [('gapB_matrix', 8, 0), ('gapB_matrix', 7, 1),('gapB_matrix', 7, 2), ('gapB_matrix', 7, 3),
                          ('align_matrix', 7, 4), ('gapB_matrix', 7, 5),('gapB_matrix', 7, 6), ('gapB_matrix', 7, 7)],
                         [('gapB_matrix', 9, 0), ('gapB_matrix', 8, 1),('gapB_matrix', 8, 2), ('gapB_matrix', 8, 3),
                          ('gapB_matrix', 8, 4), ('align_matrix', 8, 5),('gapA_matrix', 8, 6), ('gapA_matrix', 8, 7)],
                         [('gapB_matrix', 10, 0), ('gapB_matrix', 9, 1),('gapB_matrix', 9, 2), ('gapB_matrix', 9, 3),
                          ('gapB_matrix', 9, 4), ('gapB_matrix', 9, 5),('align_matrix', 9, 6), ('gapA_matrix', 9, 7)]],
                        dtype=object)    
    observed = NW2._back_B
    for i in range(0,len(seq3)):
        for j in range(0,len(seq4)):
            assert expected[i,j][0] == observed[i,j][0]
            assert expected[i,j][1] == observed[i,j][1]
            assert expected[i,j][2] == observed[i,j][2]
    
    
    
    #_back_A check
    expected = np.array([[('gapA_matrix', 0, 0), ('gapA_matrix', 0, 1),('gapA_matrix', 0, 2), ('gapA_matrix', 0, 3),
                          ('gapA_matrix', 0, 4), ('gapA_matrix', 0, 5),('gapA_matrix', 0, 6), ('gapA_matrix', 0, 7)],
                         [('gapA_matrix', 1, 0), ('gapB_matrix', 1, 0),('align_matrix', 1, 1), ('gapA_matrix', 1, 2),
                          ('gapA_matrix', 1, 3), ('gapA_matrix', 1, 4),('gapA_matrix', 1, 5), ('gapA_matrix', 1, 6)],
                         [('gapA_matrix', 2, 0), ('gapB_matrix', 2, 0),('gapB_matrix', 2, 1), ('align_matrix', 2, 2),
                          ('gapA_matrix', 2, 3), ('gapA_matrix', 2, 4),('gapA_matrix', 2, 5), ('gapA_matrix', 2, 6)],
                         [('gapA_matrix', 3, 0), ('gapB_matrix', 3, 0),('gapB_matrix', 3, 1), ('gapB_matrix', 3, 2),
                          ('align_matrix', 3, 3), ('gapA_matrix', 3, 4),('gapA_matrix', 3, 5), ('gapA_matrix', 3, 6)],
                         [('gapA_matrix', 4, 0), ('gapB_matrix', 4, 0),('gapB_matrix', 4, 1), ('align_matrix', 4, 2),
                          ('gapB_matrix', 4, 3), ('align_matrix', 4, 4),('gapA_matrix', 4, 5), ('align_matrix', 4, 6)],
                         [('gapA_matrix', 5, 0), ('gapB_matrix', 5, 0),('gapB_matrix', 5, 1), ('align_matrix', 5, 2),
                          ('gapA_matrix', 5, 3), ('gapA_matrix', 5, 4),('align_matrix', 5, 5), ('gapA_matrix', 5, 6)],
                         [('gapA_matrix', 6, 0), ('gapB_matrix', 6, 0),('gapB_matrix', 6, 1), ('gapB_matrix', 6, 2),
                          ('align_matrix', 6, 3), ('gapA_matrix', 6, 4),('gapA_matrix', 6, 5), ('align_matrix', 6, 6)],
                         [('gapA_matrix', 7, 0), ('gapB_matrix', 7, 0),('gapB_matrix', 7, 1), ('gapB_matrix', 7, 2),
                          ('align_matrix', 7, 3), ('align_matrix', 7, 4),('gapA_matrix', 7, 5), ('gapA_matrix', 7, 6)],
                         [('gapA_matrix', 8, 0), ('gapB_matrix', 8, 0),('gapB_matrix', 8, 1), ('align_matrix', 8, 2),
                          ('gapB_matrix', 8, 3), ('gapB_matrix', 8, 4),('align_matrix', 8, 5), ('gapA_matrix', 8, 6)],
                         [('gapA_matrix', 9, 0), ('gapB_matrix', 9, 0),('gapB_matrix', 9, 1), ('align_matrix', 9, 2),
                          ('gapB_matrix', 9, 3), ('gapB_matrix', 9, 4),('align_matrix', 9, 5), ('align_matrix', 9, 6)],
                         [('gapA_matrix', 10, 0), ('gapB_matrix', 10, 0),('gapB_matrix', 10, 1), ('align_matrix', 10, 2),
                          ('gapB_matrix', 10, 3), ('gapB_matrix', 10, 4),('gapB_matrix', 10, 5), ('gapB_matrix', 10, 6)]],
                        dtype=object)
    observed = NW2._back_A
    for i in range(0,len(seq3)):
        for j in range(0,len(seq4)):
            assert expected[i,j][0] == observed[i,j][0]
            assert expected[i,j][1] == observed[i,j][1]
            assert expected[i,j][2] == observed[i,j][2]
    
    




