import math
from scipy.special import logsumexp


def viterbi(emission_matrix, transmission_matrix):
    """
    :param emission_matrix:
    :param transmission_matrix:
    :return:
    """

    # Creating variables:
    # initialprob: first row of trans-mtarix
    initialprob = transmission_matrix[0]

    # Change NaN to -inf
    for i in range(len(emission_matrix)):
        if emission_matrix[i] == -1:
            continue
        for j in range(len(emission_matrix[i])):
            if emission_matrix[i][j] == "NaN":
                emission_matrix[i][j] = -float("Inf")

    delta = []
#    tilde = []   bruachst du das?

    # Run Viterbi
    # Important:  if in Ri is no list, but a == -1 -> skip this part.
    #            if skip part because of -1: use initialprob for first element != -1 Ri
    for i in range(len(emission_matrix)):

        # Case: -1 / skip
        if emission_matrix[i] == -1:
            delta.append(-1)
            continue

        # Case: Initiation:
        #       i == 0 or R[i] != -1 and R[i-1] == -1
        if i == 0 or emission_matrix[i] != -1 and emission_matrix[i - 1] == -1:
            temp_list = []
            temp = 0
            sum = 0
            for ii in range(len(emission_matrix[i])):
                temp = initialprob[ii] * math.exp(emission_matrix[i][ii])
                temp_list.append(temp)
                sum = sum + temp

            for ii in range(len(emission_matrix[i])):
                temp_list[ii] = temp_list[ii] / sum
            delta.append(temp_list)

        # Case: Consecutive:
        # Create matrix from delta-1
        # this is important!
        else:
            delta_matrix = []
            true_delta_matrix = []
            true_temp = 0
            true_temp_list = []
            for y in range(30):
                delta_matrix.append(delta[i - 1])

            # Matrix Calculation:
            temp_list = []
            temp = 0

            # Delta-Matrix .* Transition-Matrix
            for y in range(30):
                for z in range(30):
                    if delta_matrix[y][z] != "NaN":
                        true_temp = float(
                            transmission_matrix[y][z]) * float(delta_matrix[y][z])
                        #print(transmission_matrix[y][z])
                        #print(delta_matrix[y][z])
                        true_temp_list.append(true_temp)
                        # print("aus matrix:", true_temp)
                        # print("delta", delta_matrix[y][z])
                        # print("So gerechnet: ", float(delta_matrix[y][z]) * float(transmission_matrix[y][z]))

                        true_temp = 0

                    else:
                        true_temp = "NaN"
                        true_temp_list.append(true_temp)
                        true_temp = 0
                true_delta_matrix.append(true_temp_list)
                true_temp_list = []

            # Get Max-Values
            temp_max = -float("Inf")
            max_list = []
            for y in range(30):
                temp_max = max(true_delta_matrix[y])
                max_list.append(temp_max)
                temp_max = -float("Inf")

            # log
            for y in range(30):
                if max_list[y] != "NaN":
                    # print(y)
                    # print(max_list[y])
                    max_list[y] = math.log(float(max_list[y]))
                    # print(max_list[y])
                else:
                    max_list[y] = -float("inf")

            # plus emission prob
            emmi_list = []
            for y in range(30):
                if "nan" != max_list:
                    emmi_list.append(
                        float(max_list[y]) + float(emission_matrix[i][y]))
            max_list = []

            # logsumexp:
            den = logsumexp(emmi_list)


            # Exp and sub. of logsumexp
            temp_list = []
            for y in range(len(emmi_list)):
                temp_list.append(math.exp(emmi_list[y] - den))

            delta.append(temp_list)

    #Get xtilde
    xtilde = []

    for i in reversed(range(len(delta))):
        # Case: skip
        if delta[i] == -1:
            xtilde.append(-1)
            continue
        # Case: Initiation
        if i == len(delta) and delta[i] != -1 or delta[i] != -1 and delta[i + 1] == -1:
            temp = 0
            temp = delta[i].index(max(delta[i]))
            xtilde.append(temp)

        else:
            trans_row = transmission_matrix[xtilde[i + 1]]
            temp = 0
            temp_list = []
            for j in range(len(trans_row)):
                temp = 0
                temp = float(delta[i][j]) * float(trans_row[j])
                temp_list.append(temp)
            xtilde.append(temp_list.index(max(temp_list)))
            temp_list = []

    xtilde.reverse()

    return xtilde
