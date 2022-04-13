#!/usr/bin/python
__author__ = "Yunzhe Jiang"
__email__ = "yunzhe.jiang@yale.edu"
__copyright__ = "Copyright 2022"
__license__ = "GPL"
__version__ = "1.0.0"

import argparse, sys
import pandas as pd
import numpy as np

def runSW(inputFile, scoreFile, outputFile, openGap = -2, extGap = -1):
    ## Read the score file
    blosum = pd.read_csv(scoreFile, header = 0, sep = '\s+')
    ## Read the input sequence file
    with open(inputFile, 'r') as f:
        seqs = f.readlines()
    
    seq1, seq2 = seqs[0].strip(), seqs[1].strip()
    seq1, seq2 = [ ' ' + seq for seq in [seq1, seq2]]
    matrix = np.zeros((len(seq2), len(seq1)))
    max_score = 0
    
    ## Fill the scoring matrix 
    for r in range(1, len(seq2)):
        for c in range(1, len(seq1)):
            score_m = matrix[r-1, c-1] + blosum[seq2[r]][seq1[c]] ## match/mismatch (diagonal) 
            score_l = np.max(matrix[r,:c] + ((np.arange(c, 0, -1) - 1) * extGap + openGap)) ## from the left
            score_u = np.max(matrix[:r,c] + ((np.arange(r, 0, -1) - 1) * extGap + openGap)) ## from the upper
            matrix[r,c] = np.max([score_m, score_l, score_u, 0])
            max_score = matrix[r,c] if max_score <= matrix[r,c] else max_score
   
    ## Traceback
    traceback = []
    r, c = divmod(np.argmax(matrix), matrix.shape[1]) ## find the maxmimum score in the matrix
    traceback.append((r,c)) ## starting point of traceback

    while matrix[r, c] > 0 and r > 1 and c > 1: ## First row and column are all zeros
        if matrix[r-1, c-1] + blosum[seq2[r]][seq1[c]] == matrix[r, c]: ## match/mismatch (diagonal)
            r, c = r-1, c-1 
        elif (matrix[r,:c] + ((np.arange(c, 0, -1) - 1) * extGap + openGap) == matrix[r,c]).any(): ## from the left
            r, c = r, np.where(matrix[r,:c] + ((np.arange(c, 0, -1) - 1) * extGap + openGap) == matrix[r,c])[0].item()
        elif (matrix[:r,c] + ((np.arange(r, 0, -1) - 1) * extGap + openGap) == matrix[r,c]).any(): ## from the upper
            r, c = np.where(matrix[:r,c] + ((np.arange(r, 0, -1) - 1) * extGap + openGap) == matrix[r,c])[0].item(), c   
        if matrix[r, c] > 0:  ## only keep elements greater than zero in the traceback list
            traceback.append((r,c))

    ## Generate alignment
    traceback = sorted(traceback) ## from the beginning
    seq1_align, seq2_align = seq1[traceback[0][1]], seq2[traceback[0][0]]
    mid = '|' if seq1_align[-1] == seq2_align[-1] else ' '

    for ind in range(len(traceback) - 1):
        seq1_gap = True; seq2_gap = True
        cur_pos, next_pos = traceback[ind], traceback[ind + 1] ## current pos. and next pos. required for calculating gaps

        if cur_pos[0] != next_pos[0]:
            seq2_align += seq2[cur_pos[0] + 1:next_pos[0] + 1]
            seq2_gap = False
        else:  ## gaps in sequence 2, determined by sequence 1
            seq2_align += '-' * (next_pos[1] - cur_pos[1])
            mid += ' ' * (next_pos[1] - cur_pos[1])

        if cur_pos[1] != next_pos[1]:
            seq1_align += seq1[cur_pos[1] + 1:next_pos[1] + 1]
            seq1_gap = False
        else:  ## gaps in sequence 1, determined by sequence 2
            seq1_align += '-' * (next_pos[0] - cur_pos[0])
            mid += ' ' * (next_pos[0] - cur_pos[0])

        if (not seq1_gap) and (not seq2_gap):
            mid = mid + '|' if seq1_align[-1] == seq2_align[-1] else mid + ' '  ## match/mistach

    seq1_align = seq1[1:traceback[0][1]] + ''.join(['(', seq1_align, ')']) + seq1[traceback[-1][1]+1:] ## compelete sequence 1
    seq2_align = seq2[1:traceback[0][0]] + ''.join(['(', seq2_align, ')']) + seq2[traceback[-1][0]+1:] ## compelete sequence 2

    if traceback[0][0] >= traceback[0][1]: ## align to the left
        seq1_align = ' ' * (traceback[0][0] - traceback[0][1]) + seq1_align 
    else:
        seq2_align = ' ' * (traceback[0][1] - traceback[0][0]) + seq2_align  

    if len(seq2) - traceback[-1][0] >= len(seq1) - traceback[-1][1]:  ## align to the right
        seq1_align += ' ' * (len(seq2) - traceback[-1][0] - len(seq1) + traceback[-1][1])
    else:
        seq2_align += ' ' * (len(seq1) - traceback[-1][1] - len(seq2) + traceback[-1][0])
        
    mid = ' ' * (np.max(traceback[0])) + mid + ' ' * (np.max([len(seq1) - traceback[-1][1], len(seq2) - traceback[-1][0]])) ## align to the left and the right
    
    ## Write the output to the file
    scoring_df = pd.DataFrame(matrix)   ## Convert scoring matrix to a dataframe
    scoring_df.index = [''] + [char for char in seq2[1:]] ## row represents sequence2 (the first index should be None)
    scoring_df.columns = [''] + [char for char in seq1[1:]] ## column reprents sequence1 (the first column should be None)
    
    with open(outputFile, 'w') as f: ## create a new file
        f.write('-----------\n|Sequences|\n-----------\n')
        f.write('sequence1\n%s\n' % seq1[1:])
        f.write('sequence2\n%s\n' % seq2[1:])
        f.write('--------------\n|Score Matrix|\n--------------\n')

    scoring_df.to_csv(outputFile, sep = '\t', float_format = '%d', mode = 'a') ## append to the existing file

    with open(outputFile,'a') as f:  ## append to the existing file
        f.write('----------------------\n|Best Local Alignment|\n----------------------\n')
        f.write('Alignment Score:%d\n' % max_score)
        f.write('Alignment Results:\n')
        f.write('%s\n' % seq1_align)
        f.write('%s\n' % mid)
        f.write('%s\n' % seq2_align)

