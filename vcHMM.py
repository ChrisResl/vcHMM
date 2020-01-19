import math
from scipy.special import logsumexp
import pysam
import argparse
from Bio import SeqIO

import numpy as np


def get_ref_fasta(file_name):
    """
    Reads a fasta-file and return the first DNA-Sequence.
    Needs: a full file-name, e.g. : ref.fa
    :return: DNA-Sequence (as String) of a given fasta-file.
    """
    record_dict = list(SeqIO.parse(file_name, "fasta"))
    return record_dict[0].seq


def get_sam(samfile):
    """
    Reads a sam-file and returns a list with startposition, readnames,
    readsequence, [cigarstring as tuples], [queryqualities] and mapquality.
    Needs: a full file-name, e.g. : ref.fa
    :return: All sorted reads from given sam-file.
    """
    samfile = pysam.AlignmentFile(samfile, "rb")
    sam = []
    for read in samfile:
        sam.append([read.reference_start, read.query_name, read.query_sequence,
                    read.cigartuples, read.mapping_quality, read.query_qualities.tolist()])
        sam.sort()
    return sam


def get_pileup(samfile, pileupposition, new_ref_index):
    """
    Get all read bases, qualities and mapq for a given reference position.
    Needs:  samfile(with all kind of modifications from before!),
            pileupposition(Reference-Positions)

    e.g.
        R[2] -> ["A", "A", "A", "C", "C", "C"], [22, 12, 23, 20, 18, 21],  [18]

    :return: bases, qualities and mapping_qualities (mapq)
    """
    bases = []
    qualities = []
    mapping_qualities = []
    for read in samfile:
        #nem_index_position = new_ref_index.index(read[0], read[0])
        nem_index_position = int(new_ref_index.get(str(read[0])))

        if nem_index_position <= pileupposition and nem_index_position + len(read[2]) > pileupposition:
            bases.append(read[2][pileupposition - nem_index_position])
            qualities.append(read[5][pileupposition - nem_index_position])
            mapping_qualities.append(read[4])

    return bases, qualities, mapping_qualities


def create_row_transition_matrix(vector_of_pre_transition_matrix, hetrate):
    """
    Calculates a row of the transition Matrix, given a
    "pre transition matrix" and a hetrate.
    Keep in Mind:
        - Pre Transition Matrix: only 4 hidden States, because
        for predicting heterozygous polymorphism,
        one needs a more advanced Transition Matrix
        - States are 1-based, python-code is 0-based
            -> transrow(30) = row_transition_matrix[29]
     :return:  A single row for the Transition Matrix
    """

    m = 30
    row_transition_matrix = [-1 for i in range(m)]

    # transrow(30) = tprobi(4) * hetrate/32;
    row_transition_matrix[29] = vector_of_pre_transition_matrix[3] * hetrate / 32

    # transrow(1) = tprobi(1) * (1-hetrate)*(1-transrow(30));
    row_transition_matrix[0] = vector_of_pre_transition_matrix[0] * \
                               (1 - hetrate) * (1 - row_transition_matrix[29])

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
        row_transition_matrix[i] = vector_of_pre_transition_matrix[1] * \
                                   hetrate / 6 * (1 - row_transition_matrix[29])

    # transrow([11,13,14]) = (tprobi(2)/3 + tprobi(3)) * hetrate/4*(1-transrow(30));
    for i in ([10, 12, 13]):
        row_transition_matrix[i] = ((vector_of_pre_transition_matrix[1] / 3 + vector_of_pre_transition_matrix[2]) *
                                    hetrate / 4 * (1 - row_transition_matrix[29]))

    # transrow(15) = tprobi(3) * (1 - hetrate) * (1 - transrow(30));
    row_transition_matrix[14] = vector_of_pre_transition_matrix[2] * \
                                (1 - hetrate) * (1 - row_transition_matrix[29])

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
    """

    """

    # size of transition matrix: m
    m = 30
    # = [[-1 for i in range(m)] for j in range(m)]
    transition_matrix = [-1 for i in range(m)]

    ### This one: MATCH
    #  (1,:)= buildTrans(tprob(1,:), hetrate);
    transition_matrix[0] = create_row_transition_matrix(
        pre_transition_matrix[0], hetrate)

    ### This one: SNP
    # transitionmatrix(2,:)= buildTrans(tprob(2,:),hetrate);
    transition_matrix[1] = create_row_transition_matrix(
        pre_transition_matrix[1], hetrate)
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

    ### This one: Full-Delition
    # transitionmatrix(15,:) = buildTrans(tprob(3,:), hetrate);
    transition_matrix[14] = create_row_transition_matrix(
        pre_transition_matrix[2], hetrate)

    ### Part-Deletions: this are all Genotypes with ONE GAP (deletion).
    #   order is NOT the same as the paper, but same as in Matlab
    # transitionmatrix(8,:) = (transitionmatrix(2,:)+transitionmatrix(15,:))/2;
    # transitionmatrix(11,:) = transitionmatrix(8,:);
    # transitionmatrix(13,:) = transitionmatrix(8,:);
    # transitionmatrix(14,:) = transitionmatrix(8,:);

    for i in ([7, 10, 12, 13]):
        transition_matrix[i] = [-1 for i in range(m)]
        for j in range(0, len(transition_matrix[1])):
            transition_matrix[i][j] = (
                                              transition_matrix[1][j] + transition_matrix[14][j]) / 2

    ### This are: INSERTIONs
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
    transition_matrix[15] = create_row_transition_matrix(
        pre_transition_matrix[3], hetrate)
    for i in range(16, 29):
        transition_matrix[i] = transition_matrix[15]

    ### This one: Invalid State
    transition_matrix[29] = [1 / 30 for i in range(m)]

    return transition_matrix


def build_emissionmatrix(upd_sam, upd_reference, ref_index):
    """
    Create the Emission Matrix.
    :param pile_up_read:
    :param pile_up_quality:
    :param updated_reference:
    :param mapq_list:
    :return: Emission Matrix
    """
    ### Creating Genotype-Vectors / Hidden-State-Vectors
    # This is a "translation" from the MATLAB-Code, from numbers to letters.
    vector_A = [["A", "A"], ["C", "C"], ["G", "G"], ["T", "T"], ["A", "C"], ["A", "G"], ["A", "T"], ["A", "-"],
                ["C", "G"], ["C", "T"], ["C", "-"], ["G", "T"], ["G", "-"], ["T", "-"], ["-", "-"]]
    vector_C = [["C", "C"], ["A", "A"], ["G", "G"], ["T", "T"], ["A", "C"], ["C", "G"], ["C", "T"], ["C", "-"],
                ["A", "G"], ["A", "T"], ["A", "-"], ["G", "T"], ["G", "-"], ["T", "-"], ["-", "-"]]
    vector_G = [["G", "G"], ["A", "A"], ["C", "C"], ["T", "T"], ["A", "G"], ["C", "G"], ["G", "T"], ["G", "-"],
                ["A", "C"], ["A", "T"], ["A", "-"], ["C", "T"], ["C", "-"], ["T", "-"], ["-", "-"]]
    vector_T = [["T", "T"], ["A", "A"], ["C", "C"], ["G", "G"], ["A", "T"], ["C", "T"], ["G", "T"], ["T", "-"],
                ["A", "C"], ["A", "G"], ["A", "-"], ["C", "G"], ["C", "-"], ["G", "-"], ["-", "-"]]
    # keep in mind: vector_Gap = vectot_A

    ematrix = []  # Length: updated_reference * 30

    ref_index_dict = get_ref_index_dict(upd_reference, ref_index)

    # This loop is for The len of reference / Run over reference.
    # For every Reference-Base: creating emission prob. for 30 states(list)
    #
    for i in range(len(upd_reference)):

        # get pileup of reads
        # get pileup of quality
        # get mapq
        pileup, pileup_qual, mapq_list = get_pileup(upd_sam, i, ref_index_dict)

        # skip loop, if len(read-pileup) is <5 or Reference-Base is a "N"
        # append -1 as indicator for this case
        if len(pileup) < 5 or upd_reference[i] == "N":
            ematrix.append(-1)
            continue

        ## Change quality score:
        # case: all values belong to gaps, problem: gaps do not have q-scores
        # quality-score -> mapq/4
        if all(elem == pileup_qual[0] for elem in pileup_qual) and [pileup_qual[p] == "-" for p in
                                                                    range(len(pileup_qual))]:
            for y in range(len(pileup_qual)):
                pileup_qual[y] = mapq_list[y] / 4

        ## Change quality score:
        # case: somegaps are given, problem: gaps do not have q-scores!
        # missing scores: means of other scores
        elif "-" in pileup_qual:
            tempt_counter = 0
            temp_val = 0
            for y in range(len(pileup_qual)):
                if not pileup_qual[y] == "-":
                    tempt_counter = tempt_counter + 1
                    temp_val = temp_val + pileup_qual[y]

            mean = temp_val / tempt_counter
            for y in range(len(pileup_qual)):
                if pileup_qual[y] == "-":
                    pileup_qual[y] = mean

        ## For Sub-Loop of genotypes: ii
        # Creating genotype list: 1 to 30 (Python: 0 to 29)
        # Keep in Mind: Genotype = VectorPart = State
        temp_list = []
        temp_value = 0

        # case: If Reference-Base at i is a Gap
        if upd_reference[i] == "-":
            vector = vector_A  # Keep in mind, Vector A == Vector for gaps!

            # In case Ref[i] == "-":  emission prob. for states 1 to 15: NaN
            for ii in range(15):
                temp_list.append("NaN")

            # list values else: 3 if-cases fro every element in pile-up
            for ii in range(15, 30):

                ## Filter:
                # filter-case: if Vector at R[i] == Gap, Genotype == Gap and >80% of reads are gaps:
                gap_counter = 0
                for z in range(len(pileup)):
                    if pileup[z] == "-":
                        gap_counter = gap_counter + 1
                if vector[ii - 15] == ["-", "-"] and gap_counter >= len(pileup) * 0.8:
                    temp_list.append(0)
                    continue

                ## Filter
                # NaN-Checking: in vector part ii is no match with any element from pile-up-list
                nan_controll = 0
                for sub in range(len(pileup)):
                    if vector[ii - 15][0] == pileup[sub] or vector[ii - 15][1] == pileup[sub]:
                        nan_controll = 1

                if nan_controll == 0:
                    temp_list.append("NaN")
                    continue

                # Loop over Pile-Up
                ## 3 possible cases:
                temp_value = 0
                for subi in range(len(pileup)):

                    # Case1: PileUp[subi] is in both Elements of Vector Part
                    if vector[ii - 15][0] == pileup[subi] and vector[ii - 15][1] == pileup[subi]:
                        #   log10(1 - 10 ^ (-qscore[subi] / 10))
                        temp_pow = np.power(10, (-pileup_qual[subi] / 10))
                        temp = np.log10(1 - temp_pow)
                        temp_value = temp_value + temp

                    # Case2: PileUp[subi] is not in both Elements of Vector Part
                    elif pileup[subi] != vector[ii - 15][0] and pileup[subi] != vector[ii - 15][1]:
                        #   -qscore[subi] / 10 +log10(0.25)
                        temp = (-pileup_qual[subi] / 10 + np.log10(0.25))
                        temp_value = temp_value + temp

                    # Case3: PileUp[subi] is in on Element of Vector Part
                    else:
                        #   log10(0.5*(1-10^(-qscore[subi]/10)) + 0.125*10^(-qsccore[subi] / 10))
                        temp_pow = np.power(10, (-pileup_qual[subi] / 10))
                        temp = 0.5 * (1 - temp_pow) + 0.125 * temp_pow
                        temp_log = np.log10(temp)
                        temp_value = temp_value + temp_log

                temp_list.append(temp_value)
            ematrix.append(temp_list)
            temp_list = []

        # Case: If Reference-Base at i is not a Gap
        else:
            # Get Vector:
            # Keep in Mind: One does not need THIS in case of
            # a gap, because at gaps, every Vector is the same
            if upd_reference[i] == "A":
                vector = vector_A
            elif upd_reference[i] == "C":
                vector = vector_C
            elif upd_reference[i] == "G":
                vector = vector_G
            elif upd_reference[i] == "T":
                vector = vector_T
            else:
                print("Critical Error at creating emission-matrix! #2349862")

            # Sub-Loop of genotypes: ii
            for ii in range(15):
                temp_value = 0

                ## Filter:
                # skip loop, if len(read-pileup) is <5 or Reference-Base is a "N"
                # append -1 as indicator for this case
                gap_counter = 0
                for z in range(len(pileup)):
                    if pileup[z] == "-":
                        gap_counter = gap_counter + 1
                if vector[ii] == ["-", "-"] and gap_counter >= len(pileup) * 0.8:
                    temp_list.append(0)
                    continue

                ## Filter
                # NaN-Checking: in vector part ii is no match with any element from pile-up-list
                nan_controll = 0
                for sub in range(len(pileup)):
                    if vector[ii][0] == pileup[sub] or vector[ii][1] == pileup[sub]:
                        nan_controll = 1

                if nan_controll == 0:
                    temp_list.append("NaN")
                    continue

                # Loop over Pile-Up
                for subi in range(len(pileup)):

                    # Case1: PileUp[subi] is in both Elements of Vector Part
                    if vector[ii - 15][0] == pileup[subi] and vector[ii - 15][1] == pileup[subi]:
                        #   log10(1 - 10 ^ (-qscore[subi] / 10))
                        temp_pow = np.power(10, (-pileup_qual[subi] / 10))
                        temp = np.log10(1 - temp_pow)
                        temp_value = temp_value + temp

                    # Case2: PileUp[subi] is not in both Elements of Vector Part
                    elif pileup[subi] != vector[ii - 15][0] and pileup[subi] != vector[ii - 15][1]:
                        #   -qscore[subi] / 10 +log10(0.25)
                        temp = (-pileup_qual[subi] / 10 + np.log10(0.25))
                        temp_value = temp_value + temp

                    # Case3: PileUp[subi] is in on Element of Vector Part
                    else:
                        #   log10(0.5*(1-10^(-qscore[subi]/10)) + 0.125*10^(-qsccore[subi] / 10))
                        temp_pow = np.power(10, (-pileup_qual[subi] / 10))
                        temp = 0.5 * (1 - temp_pow) + 0.125 * temp_pow
                        temp_log = np.log10(temp)
                        temp_value = temp_value + temp_log

                temp_list.append(temp_value)

            # In case: Ref[i] != "-": emission prob. for states 15-30: nan
            for ii in range(15, 30):
                temp_list.append("NaN")
            ematrix.append(temp_list)
            temp_list = []

    return ematrix


def viterbi(emission_matrix, transmission_matrix):
    """
    Does Vitervi Hidden Markov for a given transition and emission matrix.
    2 Parts.
    -Viterbi
    -Xtilde

    :param emission_matrix:
    :param transmission_matrix:
    :return: xtilde (path of most likely Hidden States)
    """

    ## Creating variables:
    # initialprob: first row of trans-matrix
    initialprob = transmission_matrix[0]
    delta = []

    # Change NaN to -inf (np.NINF)
    for i in range(len(emission_matrix)):
        if emission_matrix[i] == -1:
            continue
        for j in range(len(emission_matrix[i])):
            if emission_matrix[i][j] == "NaN":
                emission_matrix[i][j] = np.NINF

    ####################
    #### Run Viterbi ###
    ####################

    # Important:  if -1 -> skip this part and.
    #             if skip because of -1: use initialprob for first next element != -1 in Ri
    for i in range(len(emission_matrix)):

        ## Case: -1 / skip
        # append -1 as indicator for -1
        # e.g.
        #   x x x -1  x x x x x x ...
        #   z z z  s  y y y y y ...
        #   z z z  k  y y y y ...
        #   z z z  i  y y y y ...
        #   z z z  p  y y y y ...
        #
        # Why: Every -1 cases cut the matrix. Each part must be calculated separately.
        if emission_matrix[i] == -1:
            delta.append(-1)
            continue

        ## Case: Initialization
        #   R[i] == 0 or (R[i] != -1 and R[i-1] == -1)
        #   This starts a separately calculation
        if i == 0 or (emission_matrix[i] != -1 and emission_matrix[i - 1] == -1):
            temp_list = []
            temp = 0
            sum = 0
            for ii in range(len(emission_matrix[i])):
                temp = initialprob[ii] * np.exp(emission_matrix[i][ii])
                temp_list.append(temp)
                sum = sum + temp

            for ii in range(len(emission_matrix[i])):
                temp_list[ii] = temp_list[ii] / sum
            delta.append(temp_list)

        # Case: Consecutive sequence
        #   -e.g.  ACGTACGT, without -1
        else:
            delta_matrix = []
            # Keep in mind: delta_matrix != delta
            # delta_matrix :         a temp. Matrix for calculation
            # true_delta_matrix :    a temp. Matrix for calculation
            # delta:                 saves solutions of needed calculations

            true_delta_matrix = []
            true_temp = 0
            true_temp_list = []
            for y in range(30):
                delta_matrix.append(delta[i - 1])

            # Matrix Calculation:
            #   - Delta-Matrix .* Transition-Matrix
            for y in range(30):
                for z in range(30):
                    if delta_matrix[y][z] != "NaN":
                        true_temp = float(
                            transmission_matrix[y][z]) * float(delta_matrix[y][z])
                        true_temp_list.append(true_temp)

                    else:
                        true_temp = "NaN"
                        true_temp_list.append(true_temp)
                true_delta_matrix.append(true_temp_list)
                true_temp_list = []

            # Get Max-Values from true_delta_matrix-rows
            # e.g.
            # [1,2,3,4] -> [4,
            # [2,3,4,5] ->  5,
            # [4,3,2,1] ->  4,
            # [5,5,5,7] ->  7]
            # Input Matrix, Output List
            temp_max = np.NINF
            max_list = []
            for y in range(30):
                temp_max = np.max(true_delta_matrix[y])
                max_list.append(temp_max)
                temp_max = np.NINF

            # Get Log-Values from Max-Values(List)
            for y in range(30):
                if max_list[y] != "NaN":
                    # print(y)
                    # print(max_list[y])
                    max_list[y] = np.log(float(max_list[y]))
                else:
                    max_list[y] = np.NINF

            # Every value(in list) plus emission prob.
            emmi_list = []
            for y in range(30):
                if "nan" != max_list[y]:  # nan to NaN
                    emmi_list.append(
                        float(max_list[y]) + float(emission_matrix[i][y]))
            max_list = []

            # Logsumexp:
            # First exp, than sum of list and log
            # get a single value
            den = logsumexp(emmi_list)

            # Exp of every value(list) and sub. of logsumexp
            temp_list = []
            for y in range(len(emmi_list)):
                temp_list.append(np.exp(emmi_list[y] - den))

            delta.append(temp_list)

    ####################
    #### Get xtilde  ###
    ####################
    # Get Path of Hidden Satets
    # Get most likely Hidden States as a List

    xtilde = []
    # reverse, bacuse xtilde is done backwards
    delta.reverse()

    for i in (range(len(delta))):
        # Case: Skip
        #   x x x -1  x x x x x x ...
        #   z z z  s  y y y y y ...
        #   z z z  k  y y y y ...
        #   z z z  i  y y y y ...
        #   z z z  p  y y y y ...
        if delta[i] == -1:
            xtilde.append(-1)
            continue

        # Case: Initialization of new Sub-Part
        if (i == 0 and delta[i] != -1) or (delta[i] != -1 and delta[i + 1] == -1):
            temp = delta[i].index(np.max(delta[i]))
            xtilde.append(temp)

        else:
            trans_row = transmission_matrix[xtilde[i - 1]]
            temp_list = []
            for j in range(len(trans_row)):
                temp = 0
                temp = float(delta[i][j]) * float(trans_row[j])
                temp_list.append(temp)
            xtilde.append(temp_list.index(np.max(temp_list)))
            temp_list = []

    xtilde.reverse()
    print("Viterbi is done.")
    return xtilde


def find_base_state(xtilde, upd_ref):
    """
    Finds the most likely genotypes / Base State at every position of updated reference.

    e.g. at R[i] Base is G and hidden state(from xtilde) is 1: -> Base State is [G, G]

    :param xtilde:  Hidden State of Genome Reference
    :param upd_ref: Genome Reference with Gaps
    :return: base_states
    """
    base_state = []

    ### Creating Genotype-Vectors / Hidden-State-Vectors
    vector_A = [["A", "A"], ["C", "C"], ["G", "G"], ["T", "T"], ["A", "C"], ["A", "G"], ["A", "T"], ["A", "-"],
                ["C", "G"], ["C", "T"], ["C", "-"], ["G", "T"], ["G", "-"], ["T", "-"], ["-", "-"]]
    vector_C = [["C", "C"], ["A", "A"], ["G", "G"], ["T", "T"], ["A", "C"], ["C", "G"], ["C", "T"], ["C", "-"],
                ["A", "G"], ["A", "T"], ["A", "-"], ["G", "T"], ["G", "-"], ["T", "-"], ["-", "-"]]
    vector_G = [["G", "G"], ["A", "A"], ["C", "C"], ["T", "T"], ["A", "G"], ["C", "G"], ["G", "T"], ["G", "-"],
                ["A", "C"], ["A", "T"], ["A", "-"], ["C", "T"], ["C", "-"], ["T", "-"], ["-", "-"]]
    vector_T = [["T", "T"], ["A", "A"], ["C", "C"], ["G", "G"], ["A", "T"], ["C", "T"], ["G", "T"], ["T", "-"],
                ["A", "C"], ["A", "G"], ["A", "-"], ["C", "G"], ["C", "-"], ["G", "-"], ["-", "-"]]

    # Keep in Mind:  Vector for Base A = Vector for Base "-"
    vector_A_ex = vector_A
    vector_A_ex.extend(vector_A)

    vector_C_ex = vector_C
    vector_C_ex.extend(vector_A)

    vector_G_ex = vector_G
    vector_G_ex.extend(vector_A)

    vector_T_ex = vector_T
    vector_T_ex.extend(vector_A)

    # Loop over updated reference:
    for i in range(len(upd_ref)):

        # Case: xtilde(i) == -1:
        #   skip this, base = -1
        if xtilde[i] == -1:
            base_state.append(-1)

        else:
            # Get Base at Ri:
            if upd_ref[i] == "A":
                base_state.append(vector_A_ex[xtilde[i]])

            elif upd_ref[i] == "C":
                base_state.append(vector_C_ex[xtilde[i]])

            elif upd_ref[i] == "G":
                base_state.append(vector_G_ex[xtilde[i]])

            elif upd_ref[i] == "T":
                base_state.append(vector_T_ex[xtilde[i]])

            elif upd_ref[i] == "N":
                base_state.append(upd_ref[i])
                print("Error. #9872536")
                # If this Error: Ns are forbidden.
                # Change logic request for selection.

            elif upd_ref[i] == "-":
                base_state.append(vector_A_ex[xtilde[i]])

            else:
                print("Critical Error. Unknown Base at Position ",
                      i, " in updated Reference! #666987")
                base_state.append(upd_ref[i])

    print("States are done.")
    return base_state


def create_variant_calling_output(ref, upd_ref, base_states, xtilde, original_index):
    variants = ""
    variant_list = []

    """
    variants = str(i) + "\t" + str(ref[i]) + "\t" + base_states[i +gap_counter][0] + "," +  base_states[i + gap_counter][1]
    """

    i = 0
    while i < len(upd_ref):

        ### Case: Skip this Element
        if xtilde[i] == 0 or xtilde[i] == 29 or xtilde[i] == -1:
            # Case: Hidden State: 1, 30 and -1
            #       1:  No Mutation
            #       30: Not valid state
            #       -1: not Data
            # -> skip this  loop-part # improve runtime
            i = i + 1
            continue

        ### Get Variables

        variant_temp = ""

        # For Difference between reference and updated reference
        # And cases ....

        # for deletions
        original_ref_position = original_index[i]
        base_at_original_ref_position = ref[original_ref_position]
        base_before_orig_ref_position = ref[original_ref_position - 1]

        # for insertions
        original_ref_position_ins = original_index[i - 1] + 1
        base_before_orig_ref_position_ins = ref[original_ref_position_ins - 1]

        original_ref_position_out = original_ref_position + 1
        before_ori_ref_position = original_ref_position

        ### Are base states equal?
        if base_states[i][0] == base_states[i][1]:
            base_states_equal = "y"
            base_state = base_states[i][0]
        else:
            base_states_equal = "n"
            base_state_1 = base_states[i][0]
            base_state_2 = base_states[i][1]

        ### Case: gap in updated reference ############################################################
        #   Insertion
        #   3 Kinds of Insertions:
        #   - one Insertion             e.g  A AG
        #   - two different Insertions  e.g. A AG, AC
        #   - multible Insertions       e.g. A AGCGCGT
        if upd_ref[i] == "-":

           ### Sub-Case 1 and 2: ##################
           # Part after this must be skiped
            if upd_ref[i + 1] != "-":

                ### Sub-Case 2: 2 different Insertions:
                #   e.g. 2345 A AG, AC
                if base_states_equal == "n":
                    if base_state_1 != "-" and base_state_2 != "-":
                        variant_temp = str(original_ref_position_ins) + "\t" + str(base_before_orig_ref_position_ins) + "\t" + str(base_before_orig_ref_position_ins) + str(base_state_1) + "," + str(base_before_orig_ref_position_ins) + str(base_state_2)
                        i = i + 1
                        variant_list.append(variant_temp)
                        continue
                    else:
                        if base_state_2 == "-":
                            not_gap = base_state_1
                        else:
                            not_gap = base_state_1
                    variant_temp = str(original_ref_position_ins) + "\t" + str(
                        base_before_orig_ref_position_ins) + "\t" + str(base_before_orig_ref_position_ins) + str(
                        not_gap)
                    i = i + 1
                    variant_list.append(variant_temp)
                    continue

                ###Sub-Case 1: 1 kind of Insertion:
                #   e.g. 2345 A AG
                else:
                    variant_temp = str(original_ref_position_ins) + "\t" + str(base_before_orig_ref_position_ins) + "\t" + str(base_before_orig_ref_position_ins) + str(base_state)
                    i = i + 1
                    variant_list.append(variant_temp)
                    continue

            ### Sub-Case 3:  #######################
            #   e.g. 2345 A ACGTCGT
            #   e.g. 2345 A ACGTCGT, ACCCCCC
            elif upd_ref[i + 1] == "-":
                ### Do a loop to get multible insertions
                loop_indicator = True
                i_start = i
                # 2 Different Base-Sequence-Elongations:
                base_seq_1 = base_states[i][0] + base_states[i + 1][0]
                base_seq_2 = base_states[i][1] + base_states[i + 1][1]
                i = i + 1
                while loop_indicator == True:
                    if upd_ref[i + 1] == "-":
                        base_seq_1 = base_seq_1 + base_states[i + 1][0]
                        base_seq_2 = base_seq_2 + base_states[i + 1][1]
                        i = i + 1
                    else:
                        loop_indicator = False

                #   e.g. 2345 A ACGTSCGT
                if base_seq_1 == base_seq_2:
                    variant_temp = str(original_ref_position_ins) + "\t" + str(base_before_orig_ref_position_ins) + "\t" + str(base_before_orig_ref_position_ins) + str(base_seq_1)
                    i = i + 1
                    variant_list.append(variant_temp)
                    continue

                #   e.g. 2345 A ACGT, ACCC
                else:
                    base_seq_1 = base_seq_1.replace("-", "")
                    base_seq_2 = base_seq_2.replace("-", "")

                    if len(base_seq_1) > 0 and len(base_seq_2) > 0:
                        variant_temp = str(original_ref_position_ins) + "\t" + str(base_before_orig_ref_position_ins) + "\t" + str(base_before_orig_ref_position_ins) + str(base_seq_1) + "," + str(base_before_orig_ref_position_ins) + str(base_seq_2)
                        i = i + 1
                        variant_list.append(variant_temp)
                        continue
                    else:
                        if len(base_seq_1) > 0:
                            temp = base_seq_1
                        else:
                            temp = base_seq_2
                        variant_temp = str(original_ref_position_ins) + "\t" + str(
                            base_before_orig_ref_position_ins) + "\t" + str(base_before_orig_ref_position_ins) + str(
                            temp)
                        i = i + 1
                        variant_list.append(variant_temp)
                        continue

                    ### Case: Not-Gap in updated Reference ###########################################################
        #   Can be:
        #       -Deletion
        #       -Part-Deletion
        #       -One SNP
        #       -Two different SNPs
        else:
            ### Case: Complete Deletion / Deletion on both Strings
            #       e.g. 2345 CG C
            if xtilde[i] == 14:
                variant_temp = str(before_ori_ref_position) + "\t" + str(base_before_orig_ref_position) + str(base_at_original_ref_position) + "\t" + str(base_before_orig_ref_position)
                i = i + 1
                variant_list.append(variant_temp)
                continue

            ### Case: Part-Deletion w/o SNP
            #   e.g. 2345 AG AT, A or 2345 AG AC, A
            #   Base_state 1 or 2 are Gaps!
            elif base_states_equal == "n" and (base_state_1 == "-" or base_state_2 == "-"):
                # get not-gap
                if base_state_2 == "-":
                    not_gap = base_state_1
                else:
                    not_gap = base_state_1
                if not_gap != base_at_original_ref_position:
                    variant_temp = str(before_ori_ref_position) + "\t" + str(base_before_orig_ref_position) + str(base_at_original_ref_position) + "\t" + str(base_before_orig_ref_position) + str(not_gap) + "," + str(base_before_orig_ref_position)
                    i = i + 1
                    variant_list.append(variant_temp)
                    continue
                else:
                    variant_temp = str(before_ori_ref_position) + "\t" + str(base_before_orig_ref_position) + str(
                        base_at_original_ref_position) + "\t" + str(base_before_orig_ref_position)
                    i = i + 1
                    variant_list.append(variant_temp)
                    continue


            ### Case: Two-SNP
            #   e.g. 2345 A C,G
            elif base_states_equal == "n":

                if base_state_1 != upd_ref[i] and base_state_2 != upd_ref[i]:
                    variant_temp = str(original_ref_position_out) + "\t" + str(base_at_original_ref_position) + "\t" + str(base_state_1) + "," + str(base_state_2)
                    i = i + 1
                    variant_list.append(variant_temp)
                    continue
                #   e.g. 2345 A G but one Base-State is equal to Reference-Base
                else:
                    # get Unequal Base State
                    if base_state_1 != upd_ref[i]:
                        temp = base_state_1
                    else:
                        temp = base_state_2
                    variant_temp = str(original_ref_position_out) + "\t" + str(base_at_original_ref_position) + "\t" + str(temp)
                    i = i + 1
                    variant_list.append(variant_temp)
                    continue


            ### Case: One-SNP
            #   e.g. 2345 A C
            elif base_states_equal == "y":
                variant_temp = str(original_ref_position_out) + "\t" + str(base_at_original_ref_position) + "\t" + str(base_state)
                i = i + 1
                variant_list.append(variant_temp)
                continue
            else:
                variant_temp = str(original_ref_position_out) + "\t" + str(base_at_original_ref_position) + "\t" + str(base_states[i][0]) + "," + str(base_states[i][1])
                print("Error. #678890123")

    return variant_list


def create_output_file(values, file):
    # add name-column to list
    new_values = []
    for element in values:
        print(element)
        temp = "simref\t" + element
        new_values.append(temp)

    # add id-column
    new_values2 = []
    id_column_value = "\t.\t"
    needle = "\t"
    for i in range(len(new_values)):
        start = new_values[i].find(needle)
        n = 2
        while start >= 0 and n > 1:
            start = new_values[i].find(needle, start + len(needle))
            n -= 1
        temp = new_values[i][0:start] + id_column_value + \
            new_values[i][start + 1:len(new_values[i])]
        new_values2.append(temp)

    # add header
    head_list = []
    meta_info = "##fileformat=VCFv4.0"
    head_list.append(meta_info)

    marker_info = "#CHROM\tPOS\tID\tREF\tALT"
    head_list.append(marker_info)

    head_list.extend(new_values2)

    # w file
    with open(file, 'w') as output:
        for element in head_list:
            output.write(str(element) + '\n')


def parser():
    """
    Allows to use the tool from the command line.
    """
    parser = argparse.ArgumentParser(
        description='Calculating Variants in fasta file')
    parser.add_argument('-i', '--input', type=str,
                        help='Input reference file', required=True)
    parser.add_argument('-r', '--reads', type=str,
                        help='Input sam read file', required=True)
    parser.add_argument('-o', '--output', type=str,
                        help='Output vcf file', required=True)

    args = parser.parse_args()

    return args


def x_read(ref, reads):

    ##############################################################################################################
    #####
    ##############################################################################################################
    upd_ref = ""
    new_reads = []

    uni_insertions, updated_reads = get_uni_insertions_and_update_reads(reads)
    upd_ref, upd_ref_index = update_ref(ref, uni_insertions)
    final_reads = update_reads(upd_ref, upd_ref_index, updated_reads, reads)

    final_reads = tuple(final_reads)

    return upd_ref, final_reads, upd_ref_index


def update_reads(upd_ref, upd_new_index, updated_reads, reads):
    """
    Insert Gaps from Insertions from other Reads!

    """
    final_reads = []

    for i in range(len(reads)):
        memory_i = ""
        #print(reads[i])
        start_pos = reads[i][0]
        read_name = reads[i][1]
        read_cigar = reads[i][3]
        mapq = reads[i][4]

        read_seq = updated_reads[i][0]
        read_qual = updated_reads[i][1]


        start_pos_new_index = upd_new_index.index(start_pos)
        memory_i = (start_pos_new_index + len(read_seq)) * 2
        ref_part = upd_ref[start_pos_new_index:memory_i]

        #
        # get cigar-like-string
        #

        cls = ""
        for element in read_cigar:
            cls = cls + str(element[0]) * element[1]

        for u in range(len(cls)):
            if (cls[u] == "0" or cls[u] == "-") and ref_part[u] == "-":
                read_seq = read_seq[0:u] + "-" + read_seq[u:len(read_seq)]
                read_qual = read_qual[0:u] + ["-"] + read_qual[u:len(read_qual)]
                cls = cls[0:u] + "-" + cls[u:len(cls)]
                #print("si")

            # if (cls[u] == "0" or cls[u] == "2") and ref_part[u] != "-":
            #     continue
            # elif cls[u] == "1" and ref_part[u] == "-":
            #     continue
            # elif cls[u] != "1" and ref_part[u] == "-" and cls[u] != "2":
            # insert gap in seq-read:
            #     print(cls[u], cls)
            #     print(ref_part[u], ref_part)
            #     read_seq = read_seq[0:u] + "-" + read_seq[u:len(read_seq)]
            #     cls = cls[0:u] + "-" + cls[u:len(cls)]
                # ref_part = ref_part + upd_ref[memory_i + u]
        temp = []
        temp.append(start_pos)
        temp.append(read_name)
        temp.append(read_seq)
        temp.append(read_cigar)
        temp.append(mapq)
        temp.append(read_qual)
        final_reads.append(temp)

    return final_reads


def update_ref(ref, uni_inse):
    ref_index = []

    for i in range(len(ref)):
        ref_index.append(i)

    for element in uni_inse:
        true_position = ref_index.index(element[0])
        # This is not enougth:
        ref = ref[0:true_position] + element[1] * "-" + ref[true_position:len(ref)]
        for i in range(element[1]):
            #ref_index = ref_index[0:true_position] + [element[0]] + ref[true_position:len(ref)]
            ref_index.insert(true_position, 0)

    ref_index = tuple(ref_index)
    ref = tuple(ref)
    return ref, ref_index


def get_uni_insertions_and_update_reads(reads):
    """
    Insert gaps in Reads, if Deletion in Cigar Strng
    Insert gaps in gap-list, if Insertion is in Cigar String

    :return:
    """
    uni_insert = []
    updated_reads = []
    for i in range(len(reads)):
        # cigar => element[3]
        #print(read)
        start_pos = reads[i][0]
        read_name = reads[i][1]
        read_cigar = reads[i][3]

        if read_cigar == None:
            del reads[i]
            continue

        mapq = reads[i][4]
        read_seq = reads[i][2]
        read_qual = reads[i][5]

        current_read_position = start_pos
        current_ref_position = start_pos

        ### Creating Cigar-Like-String:   cls
        cls = ""
        for element in read_cigar:
            cls = cls + str(element[0]) * element[1]

        ### Get Start-Value for insertion-kaskade
        kaskade_check = 0
        insertions_len = 0

        ### Checking Sequence with Cigar-Like-String
        for j in range(len(cls)):
            current_read_position = current_read_position + 1
            current_ref_position = current_ref_position + 1

            ### Case M (0):
            #   Skip this part
            if cls[j] == "0":
                if kaskade_check == 1:
                    # add locations anf length of insertions to insertions_list
                    uni_insert.append((ref_gap_location, insertions_len))
                    kaskade_check = 0
                    insertions_len = 0
            ### Case I (1): Insertions
            #
            if cls[j] == "1":
                current_ref_position = current_ref_position -1
                # get gap location:
                if kaskade_check == 0:
                    ref_gap_location = current_ref_position

                # get len of insertion
                # problem:  only terminal i-position give
                #           information about true len. of insertions
                #           for this case: if condition
                #           - Insertions Starts in cls[j] == 1
                #           - Insertions Ends in cls[j] != 1
                kaskade_check = 1
                insertions_len = insertions_len + 1

            ### Case D (2): Deletion
            if cls[j] == "2":
                #   Gap in read_seq -> Elongation of read_seq
                read_seq = read_seq[0:j] + "-" + read_seq[j:len(read_seq)]
                #   Gap in Read_qual
                read_qual = read_qual[0:j] + ["-"] + read_qual[j:len(read_qual)]

                if kaskade_check == 1:
                    # add locations anf length of insertions to insertions_list
                    uni_insert.append((ref_gap_location, insertions_len))
                    kaskade_check = 0
                    insertions_len = 0

        updated_reads.append([read_seq, read_qual])

    uni_insert_set = list(dict.fromkeys(uni_insert))
    #print(uni_insert_set)
    return uni_insert_set, updated_reads


def get_ref_index_dict(ref, ref_index):
    """
    For better runtime of Pile-Up

    :param ref:
    :param ref_index:
    :return:
    """
    a = ref_index   # e.g.  (1, 2, 3, 3, 3, 4, 5)
    b = list(ref)         # e.g.  ["A", "C", "G", "-", "-", "T", "A"]

    counter = 0
    temp = 0
    c = []

    d = []
    for i in range(len(a)):
        if i == 0:
            c.append(str(a[i]))
            d.append(i)
        else:
            if a[i] == 0 and a[i - 1] != 0:
                temp = a[i - 1]
                counter = counter + 1
                c.append(str(temp) + "," + str(counter))
                d.append(i)

            elif a[i] == 0 and a[i - 1] == 0:
                counter = counter + 1
                c.append(str(temp) + "," + str(counter))
                d.append(i)

            else:
                c.append(str(a[i]))
                d.append(i)

    # >>> dictionary = dict(zip(keys, values))
    ref_dic = dict(zip(c, d))

    return ref_dic


def main():
    #### Create Variables:
    ## T-Matrix:
    # This one is pre-transition-matrix for simulated data:
    # header:                          MATCH, SNP, Deletion,  Insertion
    pre_transition_matrix_simulated = [[0.988, 0.008, 0.002, 0.002],
                                       [0.53, 0.45, 0.01, 0.01],
                                       [0.70, 0.15, 0.15, 0.0],
                                       [0.70, 0.15, 0.0, 0.15]]

    # This one is pre-transitions-matrix for real data
    # header:                       MATCH,      SNP,        Deletion,     Insertion
    pre_transition_matrix_real = [[0.9838043, 0.01474720, 0.0006085089, 0.0008400445],
                                  [0.9499207, 0.04640025,
                                   0.0014855172, 0.0021934910],
                                  [0.2879631, 0.01089283,
                                   0.6994015911, 0.0017424552],
                                  [0.4163771, 0.01984721, 0.0040161923, 0.5597594535]]

    ## Heterozygous Rate
    # For simulated data:
    hetrate_simulated = 0.01
    # For real data
    hetrate_real = 0.001


    ### Get data
    args = parser()
    sam = get_sam(args.reads)
    # sam = get_sam('data/test_10X_Coverage/read_sort.sam')
    ref_seq = get_ref_fasta(args.input)
    # ref_seq = get_ref_fasta('data/test_10X_Coverage/ref.fa')

    upd_ref, upd_reads, upd_ref_index = x_read(ref_seq, sam)

    ### Testing Pile-Up
    #ref_index_dict = get_ref_index_dict(upd_ref, upd_ref_index)
    #for i in range(2260, 3368):
    #    bases, qualities, mapping_qualities = get_pileup(upd_reads, i, ref_index_dict)
    #    print(i+1, upd_ref[i], str(upd_ref_index[i]+1), bases)


    #### Viterbi
    trans_matrix = create_transition_matrix(
        pre_transition_matrix_simulated, hetrate_simulated)
    emission_matrix = build_emissionmatrix(upd_reads, upd_ref, upd_ref_index)
    xtilde = viterbi(emission_matrix, trans_matrix)
    hidden_states = find_base_state(xtilde, upd_ref)

    #### Varient Output
    output = create_variant_calling_output(
        ref_seq, upd_ref, hidden_states, xtilde, upd_ref_index)
    create_output_file(output, args.output)


if __name__ == '__main__':
    main()
