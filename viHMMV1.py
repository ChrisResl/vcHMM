import pysam
from Bio import SeqIO


def letters_to_number_list(sequence):
    """
    This code split string into list and replaces the Letters:
    - A C G T N
    With numbers:
    - 1 2 3 4 7

    :param sequence:
    :return:List of Number
    """
    sequence_letter_list = list(sequence)
    sequence_number_list = []

    print(sequence_letter_list)

    for i in range(0, len(sequence)):
        if sequence_letter_list[i] == "A":
            sequence_number_list.append(1)
        elif sequence_letter_list[i] == "C":
            sequence_number_list.append(2)
        elif sequence_letter_list[i] == "G":
            sequence_number_list.append(3)
        elif sequence_letter_list[i] == "T":
            sequence_number_list.append(4)
        elif sequence_letter_list[i] == "-":
            sequence_number_list.append(5)
        elif sequence_letter_list[i] == "N":
            sequence_number_list.append(7)
        else:
            print("Error. Unknown Letter in given Sequence.")
    return sequence_number_list


def get_ref_fasta(file_name):
    """
    This code reads a fasta-file and return the first sequence.
    Needs: a full file-name, e.g. : ref.fa
    :return: DNA-Sequence (as String) of a given fasta-file
    """
    record_dict = list(SeqIO.parse(file_name, "fasta"))
    return record_dict[0].seq


def get_sam():
    """
    THis code reads a sam-file and return an object with all entries of it.
    :return: All reads from given sam-file.
    """
    samfile = pysam.AlignmentFile("example.sam", "rb")
    return samfile
    # for read in samfile.fetch():
    #   print(read)


def get_column(matrix, i):
    return [row[i] for row in matrix]


def create_row_transition_matrix(vector_of_pre_transition_matrix, hetrate):
    """
    Important: this code is done, works fine and is valid. FINGER WEG.
    This code calculates a row for the transition Matrix, given a pre transition matrix and a hetrate.

    - Pre Transition Matrix: only 4 hidden States, because we are predicting heterozygous polymorphism,
     one needs a more advanced Transition Matrix

     - Hetrate: A needed value, important for calculation of the advanced Transition Matrix

     :return:  A single row for the advanced Transition Matrix
    """

    ###### WICHHTIG: Immer daran denken, indizie von MATLAB MÜSSEN in python um 1 verringert werden!
    ###### Matlab rechnet ab 1, python ab 0!!!

    # size of this: 30
    m = 30
    row_transition_matrix = [-1 for i in range(m)]

    print(vector_of_pre_transition_matrix)
    print(hetrate)

    # transrow(30) = tprobi(4) * hetrate/32;
    row_transition_matrix[29] = vector_of_pre_transition_matrix[3] * hetrate / 32

    # transrow(1) = tprobi(1) * (1-hetrate)*(1-transrow(30));
    row_transition_matrix[0] = vector_of_pre_transition_matrix[0] * (1 - hetrate) * (1 - row_transition_matrix[29])

    # transrow(2: 4) = tprobi(2) * (1 - hetrate) / 3 * (1 - transrow(30));
    for i in range(1, 3 + 1):
        row_transition_matrix[i] = vector_of_pre_transition_matrix[1] * (1 - hetrate) / 3 * (
                1 - row_transition_matrix[29])

    # transrow(5: 7) = (tprobi(1) + tprobi(2) / 3) * hetrate / 4 * (1 - transrow(30));
    for i in range(4, 6 + 1):
        row_transition_matrix[i] = (vector_of_pre_transition_matrix[0] + vector_of_pre_transition_matrix[1] / 3) * \
                                   hetrate / 4 * (1 - row_transition_matrix[29])

    # transrow(8) = (tprobi(1) + tprobi(3)) * hetrate/4*(1-transrow(30));
    row_transition_matrix[7] = (vector_of_pre_transition_matrix[0] + vector_of_pre_transition_matrix[2]) * \
                               hetrate / 4 * (1 - row_transition_matrix[29])

    # transrow([9:10,12]) = tprobi(2) * hetrate/6*(1-transrow(30));
    for i in ([8, 9, 11]):
        row_transition_matrix[i] = vector_of_pre_transition_matrix[1] * hetrate / 6 * (1 - row_transition_matrix[29])

    # transrow([11,13,14]) = (tprobi(2)/3 + tprobi(3)) * hetrate/4*(1-transrow(30));
    for i in ([10, 12, 13]):
        row_transition_matrix[i] = ((vector_of_pre_transition_matrix[1] / 3 + vector_of_pre_transition_matrix[2]) *
                                    hetrate / 4 * (1 - row_transition_matrix[29]))

    # transrow(15) = tprobi(3) * (1 - hetrate) * (1 - transrow(30));
    row_transition_matrix[14] = vector_of_pre_transition_matrix[2] * (1 - hetrate) * (1 - row_transition_matrix[29])

    # transrow(16:19) = tprobi(4) * (1-hetrate) / 4 * (1 - transrow(30));
    for i in range(15, 18 + 1):
        row_transition_matrix[i] = vector_of_pre_transition_matrix[3] * (1 - hetrate) / 4 * (1 -
                                                                                             row_transition_matrix[29])

    # transrow([20:22,24,25,27]) = tprobi(4) * hetrate/8*(1-transrow(30));
    for i in ([19, 20, 21, 23, 24, 26]):
        row_transition_matrix[i] = vector_of_pre_transition_matrix[3] * hetrate / 8 * (1 -
                                                                                       row_transition_matrix[29])

    # transrow([23,26,28,29]) = tprobi(4) * hetrate/16*(1-transrow(30));
    for i in ([22, 25, 27, 28]):
        row_transition_matrix[i] = vector_of_pre_transition_matrix[3] * hetrate / 16 * (1 -
                                                                                        row_transition_matrix[29])

    return row_transition_matrix


def create_transition_matrix(pre_transition_matrix, hetrate):
    # size of transition matrix: m
    m = 30
    transition_matrix = [-1 for i in range(m)]  # = [[-1 for i in range(m)] for j in range(m)]

    ### % This one: MATCH
    #  (1,:)= buildTrans(tprob(1,:), hetrate);
    transition_matrix[0] = create_row_transition_matrix(pre_transition_matrix[0], hetrate)

    #### % This one: SNP
    # transitionmatrix(2,:)= buildTrans(tprob(2,:),hetrate);
    transition_matrix[1] = create_row_transition_matrix(pre_transition_matrix[1], hetrate)
    # transitionmatrix(3,:) = transitionmatrix(2,:);
    # transitionmatrix(4,:) = transitionmatrix(2,:);
    # transitionmatrix(5,:) = transitionmatrix(2,:);
    # transitionmatrix(6,:) = transitionmatrix(2,:);
    # transitionmatrix(7,:) = transitionmatrix(2,:);
    # transitionmatrix(9,:) = transitionmatrix(2,:);
    # transitionmatrix(10,:) = transitionmatrix(2,:);
    # transitionmatrix(12,:) = transitionmatrix(2,:);
    for i in ([2, 3, 4, 5, 6, 8, 9, 11]):
        transition_matrix[i] = transition_matrix[1]

    ### % This one: DELETE
    # transitionmatrix(15,:) = buildTrans(tprob(3,:), hetrate);
    transition_matrix[14] = create_row_transition_matrix(pre_transition_matrix[2], hetrate)

    ### % Deletions:  % this are alle Genotypes with ONE GAP (deletion). order is NOT the same as the paper, srly why?
    # transitionmatrix(8,:) = (transitionmatrix(2,:)+transitionmatrix(15,:))/2;
    # transitionmatrix(11,:) = transitionmatrix(8,:);
    # transitionmatrix(13,:) = transitionmatrix(8,:);
    # transitionmatrix(14,:) = transitionmatrix(8,:);

    for i in ([7, 10, 12, 13]):
        transition_matrix[i] = [-1 for i in range(m)]
        for j in range(0, len(transition_matrix[1])):
            transition_matrix[i][j] = (transition_matrix[1][j] + transition_matrix[14][j]) / 2

    ### % This one: INSERTION
    # transitionmatrix(16,:) = buildTrans(tprob(4,:),hetrate);
    # transitionmatrix(17,:) = transitionmatrix(16,:);
    # transitionmatrix(18,:) = transitionmatrix(16,:);
    # transitionmatrix(19,:) = transitionmatrix(16,:);
    # transitionmatrix(20,:) = transitionmatrix(16,:);
    # transitionmatrix(21,:) = transitionmatrix(16,:);
    # transitionmatrix(22,:) = transitionmatrix(16,:);
    # transitionmatrix(23,:) = transitionmatrix(16,:);
    # transitionmatrix(24,:) = transitionmatrix(16,:);
    # transitionmatrix(25,:) = transitionmatrix(16,:);
    # transitionmatrix(26,:) = transitionmatrix(16,:);
    # transitionmatrix(27,:) = transitionmatrix(16,:);
    # transitionmatrix(28,:) = transitionmatrix(16,:);
    # transitionmatrix(29,:) = transitionmatrix(16,:);
    transition_matrix[15] = create_row_transition_matrix(pre_transition_matrix[3], hetrate)
    for i in range(16, 29):
        transition_matrix[i] = transition_matrix[15]

    ### This is invalid state
    transition_matrix[29] = [1 / 30 for i in range(m)]

    return transition_matrix


# DATA for pre-transition-matrix / transition-matrix for 4 given hidden states
# this is needed for hidden states Zi i -> {1, ... , 30}
# header: MATCH, SNP, Deletion,  Insertion

# This one is pre-transition-matrix for example
# header: MATCH, SNP, Deletion,  Insertion
pre_transition_matrix_simulated = [[0.988, 0.008, 0.002, 0.002],
                                   [0.53, 0.45, 0.01, 0.01],
                                   [0.70, 0.15, 0.15, 0.0],
                                   [0.70, 0.15, 0.0, 0.15]]

# This one is pre-transitions-matrix for real data
# header: MATCH, SNP, Deletion,  Insertion
pre_transition_matrix_real = [[0.9838043, 0.01474720, 0.0006085089, 0.0008400445],
                              [0.9499207, 0.04640025, 0.0014855172, 0.0021934910],
                              [0.2879631, 0.01089283, 0.6994015911, 0.0017424552],
                              [0.4163771, 0.01984721, 0.0040161923, 0.5597594535]]

# Heterozygous Rate   --->>> why do we need this?
# For simulated data
hetrate_simulated = 0.01
# For real data
hetrate_real = 0.001

################################# DO STUFF ###########################################################

# create_transition_matrix(pre_transition_matrix_simulated, hetrate_simulated)
# for i in ([1, 2, 3, 17]):
#    print(i)

# ref_sequence = get_ref_fasta("ref.fa")

# reads_from_sam = get_sam()

# for read in reads_from_sam.fetch():
#    print(read)


a = "ACGTN"

print(letters_to_number_list(get_ref_fasta("ref.fa")))
