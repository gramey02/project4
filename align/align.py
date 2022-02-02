# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        This function reads in two sequences of letters as strings.
        The method will compute the lowest alignment score for the strings and find the corresponding sequence 
        alignments based on initialized gap penalty and gap extension values.
        Global alignment is being performed, so the Needleman-Wunsch algorithm is implemented.
        Uses the _backtrace() method to get final score and alignments as outputs.
        
        Parameters:
            seqA: str
                sequence of letters to be aligned to sequence B
            seqB: str
                sequence of letters to be aligned to sequence A
        Returns:
            self._backtrace()
                gives a Tuple[float,str,str] of the max alignment score, seqA alignment, and seqB alignment (complete with any matches/gaps)
        """
        
        # Initialize 6 matrix private attributes for use in alignment
        # create matrices for alignment scores and gaps
        self._align_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapA_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapB_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # create matrices for pointers used in backtrace procedure
        # adding dtype=np.ndarray will allow us to store objects (like tuples), rather than just integers, in the matrices
        self._back = np.ones((len(seqA) + 1, len(seqB) + 1), dtype = object) * -np.inf
        self._back_A = np.ones((len(seqA) + 1, len(seqB) + 1), dtype = object) * -np.inf
        self._back_B = np.ones((len(seqA) + 1, len(seqB) + 1), dtype = object) * -np.inf

        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        #----------------------Global Sequence Alignment-------------------------#
        
        #initialize the first row and column of the alignment/gap matrices
        self._align_matrix[0,0] = 0 #initialize the top left entry of align_matrix to zero
        for i in range(0,len(self._seqA)+1):
            self._gapB_matrix[i,0] = self.gap_open + self.gap_extend*i #initialize first column of seqA_align
        for i in range(0,len(self._seqB)+1):
            self._gapA_matrix[0,i] = self.gap_open + self.gap_extend*i #initialize first row of seqB_align
        
        #initialize the first row and column of the backtrace matrices (excluding the top left cell)
        for i in range(0,len(self._seqA)+1):
            self._back_A[i,0] = ("gapA_matrix", i, 0)
        for i in range(0,len(self._seqB)+1):
            self._back_A[0,i] = ("gapA_matrix", 0, i)
        for i in range(0,len(self._seqA)+1):
            self._back_B[i,0] = ("gapB_matrix", i, 0)
        for i in range(0,len(self._seqB)+1):
            self._back_B[0,i] = ("gapB_matrix", 0, i)
        for i in range(0,len(self._seqA)):
            self._back[i,0] = ("align_matrix", i, 0)
        for i in range(0,len(self._seqB)):
            self._back[0,i] = ("align_matrix", 0, i)
        
        #matrices have been initialized, now fill each one in
        for i in range(1,len(self._seqB)+1):
            for j in range(1,len(self._seqA)+1):
                    
                #fill alignment/gap matrices and back matrices
                
                #fill _align_matrix and _back matrix:
                #---------------------------------------------------------------------------------------------
                self._align_matrix[j,i] = self.sub_dict[(self._seqA[j-1],self._seqB[i-1])] + max(self._align_matrix[j-1,i-1],
                                                                                                 self._gapA_matrix[j-1,i-1],
                                                                                                 self._gapB_matrix[j-1,i-1])
                
                if max(self._align_matrix[j-1,i-1], self._gapA_matrix[j-1,i-1], self._gapB_matrix[j-1,i-1]) == self._align_matrix[j-1,i-1] :
                    self._back[j,i] = ("align_matrix", j-1, i-1) #where the jth, ith value of M came from, tuple(matrix, row, col)
                    
                elif max(self._align_matrix[j-1,i-1], self._gapA_matrix[j-1,i-1], self._gapB_matrix[j-1,i-1]) == self._gapB_matrix[j-1,i-1] :
                    self._back[j,i] = ("gapB_matrix", j-1, i-1) #where the jth, ith value of M came from, tuple(matrix, row, col)                    
                    
                elif max(self._align_matrix[j-1,i-1], self._gapA_matrix[j-1,i-1], self._gapB_matrix[j-1,i-1]) == self._gapA_matrix[j-1,i-1] :
                    self._back[j,i] = ("gapA_matrix", j-1, i-1) #where the jth, ith value of M came from, tuple(matrix, row, col)
                


                #fill _gapB_matrix and _back_B matrix:
                #----------------------------------------------------------------------------------------------
                self._gapB_matrix[j,i] = max(self.gap_open + self.gap_extend + self._align_matrix[j-1, i],
                                             self.gap_extend + self._gapB_matrix[j-1,i],
                                             self.gap_open + self.gap_extend + self._gapA_matrix[j-1, i])
                
                if max(self.gap_open + self.gap_extend + self._align_matrix[j-1, i],
                       self.gap_extend + self._gapB_matrix[j-1,i],
                       self.gap_open + self.gap_extend + self._gapA_matrix[j-1, i]) == self.gap_open + self.gap_extend + self._align_matrix[j-1, i] :
                    self._back_B[j,i] = ("align_matrix", j-1, i)
                    
                elif max(self.gap_open + self.gap_extend + self._align_matrix[j-1, i],
                         self.gap_extend + self._gapB_matrix[j-1,i],
                         self.gap_open + self.gap_extend + self._gapA_matrix[j-1, i]) == self.gap_open + self.gap_extend + self._gapA_matrix[j-1, i] :
                    self._back_B[j,i] = ("gapA_matrix", j-1, i)
                
                elif max(self.gap_open + self.gap_extend + self._align_matrix[j-1, i],
                         self.gap_extend + self._gapB_matrix[j-1,i],
                         self.gap_open + self.gap_extend + self._gapA_matrix[j-1, i]) == self.gap_extend + self._gapB_matrix[j-1,i] :
                    self._back_B[j,i] = ("gapB_matrix", j-1, i)
                    

                
                
                
                #fill _gapA_matrix and _back_A matrix:
                #-------------------------------------------------------------------------------------------
                self._gapA_matrix[j,i] = max(self.gap_open + self.gap_extend + self._align_matrix[j,i-1],
                                             self.gap_open + self.gap_extend + self._gapB_matrix[j,i-1],
                                             self.gap_extend + self._gapA_matrix[j,i-1])
                #fill backtrace matrices by deterimining where the current align_matrix[j,i] value came from
                if max(self.gap_open + self.gap_extend + self._align_matrix[j,i-1],
                       self.gap_open + self.gap_extend + self._gapB_matrix[j,i-1],
                       self.gap_extend + self._gapA_matrix[j,i-1]) == self.gap_open + self.gap_extend + self._align_matrix[j,i-1] :
                    self._back_A[j,i] = ("align_matrix", j, i-1)
                    
                elif max(self.gap_open + self.gap_extend + self._align_matrix[j,i-1],
                         self.gap_open + self.gap_extend + self._gapB_matrix[j,i-1],
                         self.gap_extend + self._gapA_matrix[j,i-1]) == self.gap_extend + self._gapA_matrix[j,i-1] :
                    self._back_A[j,i] = ("gapA_matrix", j, i-1)
                    
                elif max(self.gap_open + self.gap_extend + self._align_matrix[j,i-1],
                         self.gap_open + self.gap_extend + self._gapB_matrix[j,i-1],
                         self.gap_extend + self._gapA_matrix[j,i-1]) == self.gap_open + self.gap_extend + self._gapB_matrix[j,i-1] :
                    self._back_A[j,i] = ("gapB_matrix", j, i-1)
                    

                
        return self._backtrace()
    
    
    def _backtrace(self) -> Tuple[float, str, str]:
        """
        # TODO Implement the traceback procedure method below
        based on the heuristic you implement in the align method.
        The traceback method should return a tuple of the alignment
        score, the seqA alignment and the seqB alignment respectively.
        """
        #find max alignment score by comparing the bottom right entry of _align_matrix, _gapB_matrix, and _gapA_matrix
        self.alignment_score = max(self._align_matrix[len(self._seqA), len(self._seqB)],
                                   self._gapB_matrix[len(self._seqA), len(self._seqB)],
                                   self._gapA_matrix[len(self._seqA), len(self._seqB)])
        
        #identify which matrix gives the max alignment score
        if max(self._align_matrix[len(self._seqA), len(self._seqB)],
               self._gapB_matrix[len(self._seqA), len(self._seqB)],
               self._gapA_matrix[len(self._seqA), len(self._seqB)]) == self._align_matrix[len(self._seqA), len(self._seqB)] :
            cur_matrix_type = "align_matrix"
        elif max(self._align_matrix[len(self._seqA), len(self._seqB)],
                 self._gapB_matrix[len(self._seqA), len(self._seqB)],
                 self._gapA_matrix[len(self._seqA), len(self._seqB)]) == self._gapB_matrix[len(self._seqA), len(self._seqB)] :
            cur_matrix_type = "gapB_matrix"
        elif max(self._align_matrix[len(self._seqA), len(self._seqB)],
                 self._gapB_matrix[len(self._seqA), len(self._seqB)],
                 self._gapA_matrix[len(self._seqA), len(self._seqB)]) == self._gapA_matrix[len(self._seqA), len(self._seqB)] :
            cur_matrix_type = "gapA_matrix"
        
        #create a dictionary that tells you which backtrack matrix to use given a current matrix type
        which_back_matrix = {"align_matrix": self._back, "gapB_matrix": self._back_B, "gapA_matrix": self._back_A}
        
        #x and y will allow us to index into the sequence strings
        x = len(self._seqB)
        y = len(self._seqA)
        
        #cur_row and cur_col are the location markers of the current bactrace matrix entry
        cur_col = len(self._seqB)
        cur_row = len(self._seqA)
            
            
        while x>0 or y>0:
            if cur_matrix_type=="align_matrix":
                self.seqA_align = self.seqA_align + self._seqA[y-1] #add the next letter in seqA (in reverse order)
                self.seqB_align = self.seqB_align + self._seqB[x-1] #add the next letter in seqB (in reverse order)
                #update x and y to reflect that a letter in each sequence was used
                #print("x: " + str(x))
                #print("y: " + str(y))
                x -= 1
                y -= 1
                #print("x: " + str(x))
                #print("y: " + str(y))
            elif cur_matrix_type=="gapB_matrix":
                self.seqA_align = self.seqA_align +  self._seqA[y-1] #[y-1] #add the next letter in seqA (in reverse order)
                self.seqB_align = self.seqB_align + "-" #add a gap in seqB
                #update y to reflect that a letter in seqA was used
                #print("x: " + str(x))
                #print("y: " + str(y))
                y -= 1
                #print("x: " + str(x))
                #print("y: " + str(y))
            elif cur_matrix_type=="gapA_matrix":
                self.seqA_align = self.seqA_align + "-" #add a gap in seqA
                self.seqB_align = self.seqB_align + self._seqB[x-1] #[x-1] #add the next letter in seqB (in reverse order)
                #update x to reflect that a letter in seqB was used
                #print("x: " + str(x))
                #print("y: " + str(y))
                x -= 1
                #print("x: " + str(x))
                #print("y: " + str(y))
                
            backtrace = which_back_matrix[cur_matrix_type] #get the backtrack matrix corresponding to the current matrix type
            #print(backtrace)
            
            if(x>0 or y>0):
                #update current_matrix_type by finding the matrix at the location where the backtrace matrix points
                #update the current matrix type (align_matrix, gapB_matrix, or gapA_matrix) to backtrack one step
                cur_matrix_type = (backtrace[cur_row, cur_col])[0] 
                temp_cur_row = (backtrace[cur_row, cur_col])[1] #create a temp variable to update cur_row later
                temp_cur_col = (backtrace[cur_row, cur_col])[2] #create a temp variable to update cur_col later
            
                #update row and columnn indices to reflect the location you want to go to in the next backtrack matrix
                cur_row = temp_cur_row
                cur_col = temp_cur_col
                
                #print("x: " + str(x))
                #print("y: " + str(y))
                #print("cur_row: " + str(cur_row))
                #print("cur_col: " + str(cur_col))
                #print("cur_matrix_type: " + cur_matrix_type)
                #print("seqA: " + self.seqA_align[::-1])
                #print("seqB: " + self.seqB_align[::-1])
                
            
        #reverse the alignment strings
        self.seqA_align = self.seqA_align[::-1]
        self.seqB_align = self.seqB_align[::-1] 
        
        return(self.alignment_score, self.seqA_align, self.seqB_align)   
    


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
