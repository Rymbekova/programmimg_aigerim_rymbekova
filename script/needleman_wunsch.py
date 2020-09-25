score_matrix = {'AA': 2, 'AC':-1, 'AT':-1, 'AG': 0,
                'CA':-1, 'CC': 2, 'CT': 0, 'CG':-1,
                'TA':-1, 'TC': 0, 'TT': 2, 'TG':-1,
                'GA': 0, 'GC':-1, 'GT':-1, 'GG': 2}

def scoring(seq1,seq2):
    score = 0
    d = - 2
    for i in range(len(seq1)):
        score_key = seq1[i] + seq2[i]
        if score_key in score_matrix:
        	return score_matrix[score_key]
        else:
        	return d

    print (score)

seq1 = 'GATTACACCCTACT'
seq2 = 'GTCGACGCACT'

scoring(seq1,seq2)




def needle(seq1,seq2):
 	M = len(seq1)+1
 	N = len(seq2)+1
 	d = -2
 	total = [[0]*N for x in range(M)]
 	
 	
 	for i in range(1, M):
 		total[i][0] = d * i
 		

 	for j in range(1, N):
 		total[0][j] = d * j
 		


 	for i in range(1, M):
 		for j in range(1, N):
 			match = total[i-1][j-1] + scoring(seq1[i-1], seq2[j-1])
 			delete= total[i-1][j] + d
 			insert = total[i][j-1] + d
 			total[i][j] = max(match, delete, insert)
 			
 			

 	align1 = ''
	align2 = ''
	i = M-1
	j = N-1

	while i > 0 and j > 0:
		curr = total[i][j]
		diag = total[i-1][j-1]
		up = total[i][j-1]
		left = total[i-1][j]

		if curr == diag + scoring(seq1[i-1], seq2[j-1]):
			align1 = align1 + seq1[i-1]
			align2 = align2 + seq2[j-1]
			i = i - 1
			j = j - 1			
			
		elif curr == up + d:
			align1 = align1 + '-'
			align2 = align2 + seq2[j-1]
			j = j-1

		elif curr == left + d:
			align1 = align1 + seq1[i-1]
			align2 = align2 + '-'
			i = i - 1
	
	while i > 0:
		align1 = align1 + seq1[i-1]
		align2 = align2 + '-'
		i = i - 1
	while j > 0:
		align1 = align1 + '-'
		align2 = align2 + seq2[j-1]
		j = j - 1
	align1 = align1[::-1]
	align2 = align2[::-1]
	return align1,align2

print needle (seq1, seq2)




