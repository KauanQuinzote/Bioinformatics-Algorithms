# Needleman-Wunsch Algorithm for Global Alignment

# Scoring scheme
match_award = 1
mismatch_penalty = -1
gap_penalty = -1
score_matrix = []

def create_score_matrix(rows, cols):
    """Creates a score matrix with given dimensions."""
    return [[0 for col in range(cols)] for row in range(rows)]

def score(a, b):
    """Returns the score based on match/mismatch."""
    if a == b:
        return match_award
    else:
        return mismatch_penalty

def needleman_wunsch(seq1, seq2):
    global score_matrix
    """Performs global alignment using Needleman-Wunsch algorithm."""
    # Create a matrix with dimensions (len(seq1)+1) x (len(seq2)+1)
    score_matrix = create_score_matrix(len(seq1) + 1, len(seq2) + 1)
    
    # Initialize the scoring matrix
    for i in range(1, len(seq1) + 1):
        score_matrix[i][0] = gap_penalty * i
    for j in range(1, len(seq2) + 1):
        score_matrix[0][j] = gap_penalty * j

    # Fill the scoring matrix
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            match = score_matrix[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1])
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)

    # Traceback to find the best alignment
    align1, align2 = '', ''
    i, j = len(seq1), len(seq2)

    while i > 0 and j > 0:
        score_current = score_matrix[i][j]
        score_diag = score_matrix[i - 1][j - 1]
        score_up = score_matrix[i][j - 1]
        score_left = score_matrix[i - 1][j]

        if score_current == score_diag + score(seq1[i - 1], seq2[j - 1]):
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1 += seq1[i - 1]
            align2 += '-'
            i -= 1
        elif score_current == score_up + gap_penalty:
            align1 += '-'
            align2 += seq2[j - 1]
            j -= 1

    while i > 0:
        align1 += seq1[i - 1]
        align2 += '-'
        i -= 1
    while j > 0:
        align1 += '-'
        align2 += seq2[j - 1]
        j -= 1

    # Reverse the sequences since they were built backwards
    align1 = align1[::-1]
    align2 = align2[::-1]

    return align1, align2

# Test the algorithm
lst = []

with open ('Alignment Result.txt', 'w') as f:
    alignment = needleman_wunsch('AGCT', 'ATGCT')
    f.write('Aligned Sequences:\n{0}\n{1}\n'.format(alignment[0], alignment[1]))
    
    for row in score_matrix:
        f.write(f'{row}\n')
    f.write("\n")

    alignment = needleman_wunsch('ACTGTC', 'ACGTC')
    f.write('Aligned sequences ACTGTC and ACGTC:\n{0}\n{1}\n'.format(alignment[0], alignment[1]))
    
    for row in score_matrix:
        f.write(f'{row}\n')
    f.write("\n")