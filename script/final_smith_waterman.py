#Smith-Waterman algorithm to find the best local alignment or alignments of two sequences
#according to the given scoring substitution matrix and a gap value. First, we build a scoring 
#matrix F with dimensions equal to lengths of each sequence +1. The first row and column values are
#set to 0, also all negative scored matrix cells are set to 0, so to make positively scored local
#alignment cells clear and visible. When filling the matrix, we align one sequence with another 
#element by element and score considering values of matching elements or adding gaps. If none is 
#positive, then the score of the cell is 0 meaning no similarities. Then we initialize the 
#traceback matrix P of the same dimensions as F with associated cells. But in P each cell is 
#an empty list we can append directions of obraining the maximum score to (can be one or several 
#directions). The highest score in matrix F and its index is stored, if several they are appended
#to the list of maximum score positions. To traceback, we start from the maximum score until it 
#reaches 0. We can obtain several best local alignments, if we have several best scores throughout
#the matrix F and since several elements inside each cell can give the highest score. In the end, 
#we generated parts of two sequences with the highest similarity according to the scoring matrix.


#build a scoring matrix F
#len(F) = len(sequence1) + 1
#len(F[0]) = len(sequence2) + 1
#build associated matrix P with the same dimensions but each cell is an empty list
#each cell F[i][j] has to maximum value between 0, matching two elements or adding gaps
#if F[i][j] = 0, then P[i][j].append 'X' so traceback stops
#if F[i][j] = match, then P[i][j].append 'diagonal' to know the direction of the obtained value
#if F[i][j] = delete, then P[i][j].append 'row'
#if F[i][j] = insert, then P[i][j].append 'column'
#start tracing back from maximum F[i][j], if several cells have maximum values append to a list
#start from each position from the list and follow according to the directions from P[i][j]
#iterate adding characters or gaps until F[i][j] = 0
#reverse strings as we are going backwards
#append to a list of other strings and return all possible best local alignments

from exam import BLOSUM52 #import sequences and the scoring matrix from a file

seq1 = "AAAAACCCAAAAAAAAAAAAAAATTT" #testing sequences
seq2 = "TTT"

def scoring(seq1,seq2): #if gaps or other characters are given
    gap = - 2
    for i in range(len(seq1)):
        score_key = seq1[i] + seq2[i]
        if score_key in BLOSUM52:  #if a pair given in the scoring matrix
            return BLOSUM52[score_key] #if yes, return the score according to BLOSUM52
        else:
            return -2 #if not, then punish with a gap



def matrix(seq1, seq2, BLOSUM52, gap):#a function to build matrix F based on SW and traceback matrix P
    M = len(seq1)+1 #number of columns
    N = len(seq2)+1 #number of rows 

    F = [[0 for j in range(N)] for i in range(M)] #initialize matrix F, where each cell is 0
    P = [[[] for j in range(N)] for i in range(M)] #initialize matrix P, where each cell is an empty list to append directions if case maximum score is achieved from more then one direction
    for i in range(1, M): #initialize the first row
        F[i][0] = 0 #set the first row values to 0

    for j in range(1, N): #initialize the first column
        F[0][j] = 0 #set the first column values to 0



    for i in range(1, M): #filling the matrix F
        for j in range(1, N):
            match = F[i - 1][j - 1] + scoring(seq1[i - 1], seq2[j - 1]) #a value of matching according to the scoring matrix
            delete = F[i - 1][j] + gap #a value of inserting a gap into the second seq
            insert = F[i][j - 1] + gap #a value of inserting a gap into the first seq

            F[i][j] = max(0, match, delete, insert) #for each cell the maximum value is chosen, if none is positive, then max is zero
            if F[i][j] == 0: P[i][j].append('X') #associate values in F with values in P, if score in F equals 0, then append 'X'
            if F[i][j] == match: P[i][j].append('D')#if it's a match, then append 'diagonal'
            if F[i][j] == delete: P[i][j].append('R')#if it's a delete, then append 'row'
            if F[i][j] == insert: P[i][j].append('C')#if it's an insert, then append 'column'

    return F, P #return two matrices


def max_indexes(F): #a function that iterates through the matrix F to find maximum value or values and stores their score and indexes
    i = len(F) - 1 #set variables as lengths of sequences
    j = len(F[0]) - 1
    best_i = i #set the best values 
    best_j = j
    best = F[i][j] #set the best score
    best_indexes = [] #set an empty list for indexes of max scores
    for i in range(len(F)):
        for j in range(len(F[0])): #iterate through F to find max
            curr = F[i][j]
            if curr > best: #if current cell score is larger than the best
                best = curr #then it becomes the maximum
                best_i = i
                best_j = j
                best_indexes = [[i,j]] #store the position of the maximum score in the list
            elif curr == best: #if several cells have maximum score
                best_indexes.append([i,j]) #also append them into the list
    
           
    return best_indexes, best #return the list of positions and scores of maximum values


def smith(seq1, seq2, P, best_indexes): #a function based on SW to find the best local alignment starting from the maximum scores found previously
    for best_index in best_indexes: #for each position having maximum score
        i,j = best_index
        template = '' #set the first seq as a template
        target = '' #set the second seq as a target
        template_list = [] #set an empty list of local alignments from the first seq
        target_list = [] #set an empty list of local alignments from the second seq

        while i > 0 and j > 0:  #until it reaches 0, we traceback from the maximum score indexes
            if 'D' in P[i][j]: #if 'diagonal' found in the list of directions from the current max score cell
                template = template + seq1[i - 1] #then add a character to both seqs
                target = target + seq2[j - 1]
                i -= 1 #decrease length in both
                j -= 1
                #P[i+1][j+1].remove('D') #remove this path from the previous cell 'marking as visited'
            

            elif 'R' in P[i][j]: #if 'row' found in the list of directions from the current max score cell
                template = template + seq1[i - 1] #then add a character to the template 
                target = target + '-' #and insert a gap into the target
                i -= 1 #decrease
                #P[i+1][j].remove('R') #remove this path from the previous cell 'marking as visited'
            


            elif 'C' in P[i][j]: #if 'column' found in the list of directions from the current max score cell
                template = template + "-" #insert a gap into the template
                target = target + seq2[j - 1] #then add a character to the target
                j -= 1 #decrease
                #P[i][j+1].remove('C') #remove this path from the previous cell 'marking as visited'

            else:
                break #if 0 is found or all paths are discovered then stop

        template = template[::-1] #reverse both
        target = target[::-1]

        if (template not in template_list) or (target not in target_list): #append different local alignments
            template_list.append(template)
            target_list.append(target)

        print (template_list)
        print (target_list)


def main():
    gap = -2
    F, P = matrix(seq1, seq2, BLOSUM52, gap)
    best_indexes,best = max_indexes(F)
    smith(seq1, seq2, P, best_indexes)
    print'Maximum score:', best, 'at', best_indexes
    

main()