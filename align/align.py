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
        # TODO: Fill in the Needleman-Wunsch Algorithm below
        to perform global sequence alignment of seqA and seqB
        and return a tuple with the following format
        (alignment score, seqA alignment, seqB alignment)
        Also, write up a docstring for this function using the
        _read_sub_matrix as an example.
        Don't forget to comment your code!
        """
        # Initialize 6 matrix private attributes for use in alignment
        # create matrices for alignment scores and gaps
        self._align_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapA_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapB_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # create matrices for pointers used in backtrace procedure
        self._back = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_A = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_B = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        # TODO Implement the global sequence alignment here-----------------
        
        #initialize the alignment/gap matrices
        self._align_matrix[0,0] = 0 #align[0,0] initialize the first entry of align_matrix to zero
        for i in range(0,len(self._seqA)+1):
            self._gapA_matrix[i,0] = self.gap_open + self.gap_extend*i #initialize first column of seqA_align
        for i in range(0,len(self._seqB)+1):
            self._gapB_matrix[0,i] = self.gap_open + self.gap_extend*i #initialize first row of seqB_align
        
        
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
                elif max(self._align_matrix[j-1,i-1], self._gapA_matrix[j-1,i-1], self._gapB_matrix[j-1,i-1]) == self._gapA_matrix[j-1,i-1] :
                    self._back[j,i] = ("gapA_matrix", j-1, i-1) #where the jth, ith value of M came from, tuple(matrix, row, col)
                elif max(self._align_matrix[j-1,i-1], self._gapA_matrix[j-1,i-1], self._gapB_matrix[j-1,i-1]) == self._gapB_matrix[j-1,i-1] :
                    self._back[j,i] = ("gapB_matrix", j-1, i-1) #where the jth, ith value of M came from, tuple(matrix, row, col)
                


                #fill _gapB_matrix and _back_B matrix:
                #----------------------------------------------------------------------------------------------
                self._gapB_matrix[j,i] = max(self.gap_start + self.gap_extend + self.align_matrix[j-1, i],
                                             self.gap_extend + self._gapB_matrix[j-1,i],
                                             self.gap_start + self.gap_extend + self._gapA_matrix[j-1, i])
                
                if max(self.gap_start + self.gap_extend + self.align_matrix[j-1, i],
                       self.gap_extend + self._gapB_matrix[j-1,i],
                       self.gap_start + self.gap_extend + self._gapA_matrix[j-1, i]) == self.gap_start + self.gap_extend + self.align_matrix[j-1, i] :
                    self._back_B[j,i] = ("align_matrix", j-1, i)
                
                elif max(self.gap_start + self.gap_extend + self.align_matrix[j-1, i],
                         self.gap_extend + self._gapB_matrix[j-1,i],
                         self.gap_start + self.gap_extend + self._gapA_matrix[j-1, i]) == self.gap_extend + self._gapB_matrix[j-1,i] :
                    self._back_B[j,i] = ("gapB_matrix", j-1, i)
                    
                elif max(self.gap_start + self.gap_extend + self.align_matrix[j-1, i],
                         self.gap_extend + self._gapB_matrix[j-1,i],
                         self.gap_start + self.gap_extend + self._gapA_matrix[j-1, i]) == self.gap_start + self.gap_extend + self._gapA_matrix[j-1, i] :
                    self._back_B[j,i] = ("gapA_matrix", j-1, i)
                
                
                
                #fill _gapA_matrix and _back_A matrix:
                #-------------------------------------------------------------------------------------------
                self._gapA_matrix[j,i] = max(self.gap_start + self.gap_extend + self._align_matrix[j,i-1],
                                             self.gap_start + self.gap_extend + self._gapB_matrix[j,i-1],
                                             self.gap_extend + self._gapA_matrix[j,i-1])
                #fill backtrace matrices by deterimining where the current align_matrix[j,i] value came from
                if max(self.gap_start + self.gap_extend + self._align_matrix[j,i-1],
                       self.gap_start + self.gap_extend + self._gapB_matrix[j,i-1],
                       self.gap_extend + self._gapA_matrix[j,i-1]) == self.gap_start + self.gap_extend + self._align_matrix[j,i-1] :
                    self._back_A[j,i] = ("align_matrix", j, i-1)
                    
                elif max(self.gap_start + self.gap_extend + self._align_matrix[j,i-1],
                         self.gap_start + self.gap_extend + self._gapB_matrix[j,i-1],
                         self.gap_extend + self._gapA_matrix[j,i-1]) == self.gap_start + self.gap_extend + self._gapB_matrix[j,i-1] :
                    self._back_A[j,i] = ("gapB_matrix", j, i-1)
                    
                elif max(self.gap_start + self.gap_extend + self._align_matrix[j,i-1],
                         self.gap_start + self.gap_extend + self._gapB_matrix[j,i-1],
                         self.gap_extend + self._gapA_matrix[j,i-1]) == self.gap_extend + self._gapA_matrix[j,i-1] :
                    self._back_A[j,i] = ("gapA_matrix", j, i-1)
                
                

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
            cur_matrix = self._align_matrix
            cur_backMatrix_type = "align_back"
            cur_backMatrix = self._back
            
        elif max(self._align_matrix[len(self._seqA), len(self._seqB)],
                 self._gapB_matrix[len(self._seqA), len(self._seqB)],
                 self._gapA_matrix[len(self._seqA), len(self._seqB)]) == self._gapB_matrix[len(self._seqA), len(self._seqB)] :
            cur_matrix_type = "gapB_matrix"
            cur_matrix = self._gapB_matrix
            cur_backMatrix_type = "gapB_back"
            cur_backMatrix = self._back_B
            
        elif max(self._align_matrix[len(self._seqA), len(self._seqB)],
                 self._gapB_matrix[len(self._seqA), len(self._seqB)],
                 self._gapA_matrix[len(self._seqA), len(self._seqB)]) == self._gapA_matrix[len(self._seqA), len(self._seqB)] :
            cur_matrix_type = "gapA_matrix"
            cur_matrix = self._gapA_matrix
            cur_backMatrix_type = "gapA_back"
            cur_backMatrix = self._back_A
            
        x = len(self._seqB)
        y = len(self._seqA)
        
        cur_col = len(self._seqB)
        cur_row = len(self._seqA)
            
            
        
        if cur_backMatrix_type == "align_back":
            self.seqA_align = self.seqA_align + self._seqA[y-1] #add the next letter in seqA
            self.seqB_align = self.seqB_align + self._seqB[x-1] #add the next letter in seqB
            x -= 1
            y -= 1
        elif cur_backMatrix_type == "gapB_back":
            self.seqA_align = self.seqA_align + self._seqA[y-1] #add the next letter in seqA
            self.seqB_align = self.seqB_align + "-" #add a gap in seqB
            y -= 1
        elif cur_backMatrix_type == "gapA_back":
            self.seqA_align = self.seqA_align + "-" #add a gap in seqA
            self.seqB_align = self.seqB_align + self._seqB[x-1]
            x -= 1
            
        # go into current back matrix and get the current alignment/gap matrix type to see what
        # the next alignment/gap in each sequence should look like
        cur_matrix_type = (cur_backMatrix[cur_row,cur_col])[0]
        temp_cur_row = (cur_backMatrix[cur_row,cur_col])[1]
        temp_cur_col = (cur_backMatrix[cur_row,cur_col])[2]
        cur_row = temp_cur_row
        cur_col = temp_cur_col
        
        if cur_matrix_type = "align_matrix":
            cur_backMatrix_type = "align_back"
        elif cur_matrix_type = "gapB_matrix":
            cur_matrix_type = "gapB_back"
        elif cur_matrix_type = "gapA_matrix":
            cur_matrix_type = "gapA_back"
        
        
        #now update 
        
        
        
        return(self.alignment_score, self.seqA_align, self.seqB_align)
    
    
    
    
    
    
    
    
    
    
    
    
    

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        # TODO Implement the traceback procedure method below
        based on the heuristic you implement in the align method.
        The traceback method should return a tuple of the alignment
        score, the seqA alignment and the seqB alignment respectively.
        """

        
        
        
        
        
        #this is an old version
        #get max alignment score from the "bottom right corner" of each of the matrices
        self.alignment_score = max(self._align_matrix[len(self._seqA)+1, len(self._seqB)+1],
                                   self._gapA_matrix[len(self._seqA)+1, len(self._seqB)+1],
                                   self._gapB_matrix[len(self._seqA)+1, len(self._seqB)+1])
        
        if max(self._align_matrix[len(self._seqA)+1, len(self._seqB)+1],
               self._gapA_matrix[len(self._seqA)+1, len(self._seqB)+1],
               self._gapB_matrix[len(self._seqA)+1, len(self._seqB)+1]) == self._align_matrix[len(self._seqA)+1, len(self._seqB)+1] :
            #if the max alignment score is given by the alignment matrix, make note of this in a tuple
            final_position = "align_matrix"
            self.seqA_align = self.seqA_align + self._seqA[len(self._seqA)-1]
            self.seqB_align = self.seqB_align + self._seqB[len(self._seqB)-1]
            #cur_tuple = ("align_matrix", len(self._seqA)-1, len(self._seqB)-1)
        elif max(self._align_matrix[len(self._seqA)+1, len(self._seqB)+1],
                 self._gapA_matrix[len(self._seqA)+1, len(self._seqB)+1],
                 self._gapB_matrix[len(self._seqA)+1, len(self._seqB)+1]) == self._gapA_matrix[len(self._seqA)+1, len(self._seqB)+1] :
            #if the max alignment score is given by the alignment matrix, make note of this in a tuple
            final_position = "gapA_matrix"
        elif max(self._align_matrix[len(self._seqA)+1, len(self._seqB)+1],
                 self._gapA_matrix[len(self._seqA)+1, len(self._seqB)+1],
                 self._gapB_matrix[len(self._seqA)+1, len(self._seqB)+1]) == self._gapB_matrix[len(self._seqA)+1, len(self._seqB)+1] :
            final_position = "gapB_matrix"
            
            
        #get alignment list    
            
        last_matrix = final_position
            
        
        x = len(self._seqA)
        y = len(self._seqB)
        
        #while there are still letters not aligned in each sequence...
        while x>0 or y>0:
            #if the last matrix was the align_matrix, then go to the _back matrix
            if last_matrix == "align_matrix":
                #add the final letter of each sequence to their respective alignment strings
                self.seqA_align = self.seqA_align + self._seqA[x-1]
                self.seqB_align = self.seqB_align + self._seqB[y-1]
                #subtract one from the length of each sequence since one letter from each seq will be used in the alignment strings
                x -= 1
                y -= 1
                #update last matrix
                last_matrix == (self._back[x][y])[0]
            #if the final matrix was the gapA_matrix, then go to the _back_A matrix
            elif last_matrix == "gapA_matrix":
                #add the final letter of seqB to its alignment string, add a '-' to the seqA alignment string
                self.seqA_align = self.seqA_align + "-"
                self.seqB_align = self.seqB_align + self._seqB[y-1]
                #subtract one from the length of seqB since one letter from it was just used in the alignment string
                y -= 1
                #update last matrix
                last_matrix = (self._back_A[x][y])[0]                
            #if the final matrix was the gapB_matrix, then go to the _back_B matrix
            elif last_matrix == "gapB_matrix":
                #add the final letter of seqA to its alignment string, add a '-' to the seqB alignment string
                self.seqA_align = self.seqA_align + self._seqA[x-1]
                self.seqB_align = self.seqB_align + "-"
                #subtract one from the length of seqA since one letter from it was just used in the alignment string
                x -= 1
                #update last matrix
                last_matrix = (self._back_B[x][y])[0]                
            
            
            
            
            
            
            
            
            
            
            
            
        #old version
        #initialize current variables
        cur_matrix = None
        cur_matrix_name = None
        cur_row = None
        cur_col = None
        alignment_list = [] #this will hold backtrace route
        
        #set current_matrix variable
        if final_position=="align_matrix":
            cur_matrix = self._align_matrix
            cur_matrix_name = "align_matrix"
        elif final_position=="gapA_matrix":
            cur_matrix = self._gapA_matrix
            cur_matrix_name = "gapA_matrix"
        elif final_position=="gapB_matrix":
            cur_matrix = self._gapB_matrix
            cur_matrix_name = "gapB_matrix"
        #set current column and row variables
        cur_row = len(self._seqA)
        cur_col = len(self._seqB)
        
        
        
        #get initial backtrace entry
        came_from = cur_matrix[cur_row, cur_col] #came_from is in the form: (matrix name, row #, col #)
        if came_from[0]=="align_matrix":
            back_matrix = self._back
        elif came_from[0]=="gapA_matrix":
            back_matrix = self._back_A
        elif came_from[0]=="gapB_matrix":
            back_matrix = self._back_B
        back_row = came_from[1]
        back_col = came_from[2]
        alignment_list.append((cur_matrix_name, cur_row, cur_col))
        
        while np.isinf(back_matrix[back_row, back_col]) == False:
            came_from = back_matrix[back_row, back_col] #get the tuple that represents the matrix, row #, and column # where the last alignment came from
            alignment_list.append(came_from) #add this latest alignment tuple to the alignment list, which will contain the backtrace route
            
            #now update the current backtrace matrix based on where the last alignment came from
            if came_from[0]=="align_matrix":
                back_matrix = self._back
            elif came_from[0]=="gapA_matrix":
                back_matrix = self._back_A
            elif came_from[0]=="gapB_matrix":
                back_matrix = self._back_B
            #update the row and column numbers for where the last alignment came from
            back_row = came_from[1]
            back_col = came_from[2]
        
        
        #at this point, alignment_list should be filled with tuples of the "order of traversal", or the backtrace route, through the alignment matrices
        
        seqA_align = ""
        seqB_align = ""
        #now construct two strings that represent the two sequences aligned
        for i in range(0,len(alignment_list)):
            cur_alignment = alignment_list[i]
            #see if, based on the matrix name stored in the ith list entry, there is a match, or a gap in one of the sequences
            if cur_alignment[0] == "align_matrix":
                seqA_align = seqA_align + self._seqA[cur_alignment[2]-1]
                seqB_align = seqB_align + self._seqB[cur_alignment[1]-1]
            elif cur_alignment[0] == "gapA_matrix":
                seqA_align = seqA_align + self._seqA[cur_alignment[2]]
            
            
            
            
            
            
            
            
            
            
        
        
        
        
        
        
        #delete the below code
        #construct the alignment strings from the last letter(or lack of letter) to the first letter(or lack of letter) - will reverse the strings at end
        while np.isinf(cur_matrix[cur_row, cur_col]) == False:
            alignment_list.append[(cur_matrix_name, cur_matrix, cur_row, cur_col)] #this list stores the backtrace order
            #now find where the current matrix, row value, and column value came from, and store those values
            if cur_matrix_name=="align_matrix":
                came_from = self._back[cur_row, cur_col]
            elif cur_matrix_name=="gapA_matrix":
                came_from = self._back_A[cur_row, cur_col]
            elif cur_matrix_name=="gapB_matrix":
                came_from = self._back_B[cur_row, cur_col]
            #now you should have a tuple with the matrix name, matrix, row #, and column # that the last entry came from
            
            
        
        # Implement this method based upon the heuristic chosen in the align method above.
        
        #this method has no input
        
        #so this method will make use of the back, back_A, and back_B matrices to trace back what the original alignment is, and you will concatenate together what the final alignments look like, along with the final score
        
        #don't forget to reverse the strings
        
        #return tuple of alignment score, seqA alignment, and seqB alignment
        return (self.alignment_score, self.seqA_align, self.seqB_align)


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
