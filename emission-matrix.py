import math


def build_emissionmatrix(upd_sam, upd_reference):
    """
    This code create the emissionsmatrix.
    :param pile_up_read:
    :param pile_up_quality:
    :param updated_reference:
    :param mapq_list:
    :return:
    """
    # This is a "translation" from the MATLAB-Code, from numbers to letters.
    vector_A = [["A", "A"], ["C", "C"], ["G", "G"], ["T", "T"], ["A", "C"], ["A", "G"], ["A", "T"], ["A", "-"],
                ["C", "G"], ["C", "T"], ["C", "-"], ["G", "T"], ["G", "-"], ["T", "-"], ["-", "-"]]
    vector_C = [["C", "C"], ["A", "A"], ["G", "G"], ["T", "T"], ["A", "C"], ["C", "G"], ["C", "T"], ["C", "-"],
                ["A", "G"], ["A", "T"], ["A", "-"], ["G", "T"], ["G", "-"], ["T", "-"], ["-", "-"]]
    vector_G = [["G", "G"], ["A", "A"], ["C", "C"], ["T", "T"], ["A", "G"], ["C", "G"], ["G", "T"], ["G", "-"],
                ["A", "C"], ["A", "T"], ["A", "-"], ["C", "T"], ["C", "-"], ["T", "-"], ["-", "-"]]
    vector_T = [["T", "T"], ["A", "A"], ["C", "C"], ["G", "G"], ["A", "T"], ["C", "T"], ["G", "T"], ["T", "-"],
                ["A", "C"], ["A", "G"], ["A", "-"], ["C", "G"], ["C", "-"], ["G", "-"], ["-", "-"]]

    ematrix = []  # Length: updated_reference * 30

    # test-values. clc this after get-read-list is ready.
    pileup_testmatrix = ["A", "A", "A", "A", "A"], ["C", "C", "C", "C", "A"], [
        "A", "A", "A", "A", "A"], ["T", "T", "T", "T", "T"], ["A", "A", "A", "A"]
    pileup_qual_testmatrix = [1, 2, 3, 4, 5, 6, 7, 8, "-"]

    mapq_list = [1, 1, 2, 3, 4, 5]

    # This loop is for The len of reference. / Run over reference.
    for i in range(len(upd_reference)):
        # get pileup of reads
        # get pileup of  quality
        # get mapq
        pileup = pileup_testmatrix[i]
        pileup_qual = pileup_qual_testmatrix
        # brauchst du hier nicht die mapq???

        # Control:
        # skip sub-loop, if read-pileup is <5 or Reference-Base is a "N"!
        if len(pileup) < 5 or upd_reference[i] == "N":
            ematrix.append(-1)
            continue

        # Change quality score:
        # case: all values belong to gaps: mapq/4
        if all(elem == pileup_qual[0] for elem in pileup_qual) and pileup_qual == "-":
            for y in range(len(pileup_qual)):
                pileup_qual[y] = mapq_list[y] / 4

        # case: gaps are given, problem: gaps do not have q-scores!
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

        # Sub-Loop for genotypes:
        # Creating genotype list: 1 to 30 (0 to 29)
        temp_list = []
        temp_value = 0

        # case: if reference-base at i is a gap
        if upd_reference[i] == "-":
            vector = vector_A  # Keep in mind, Vector A == Vector for gaps!

            # in case gap: genotype values 1 to 15: NaN
            for ii in range(15):
                temp_list.append("NaN")

            # list values else: 3 if-cases fro every element in pile-up
            for ii in range(15, 30):

                # Checking: case: if Vector at Gap, Gap and >80% of reads are gaps:
                gap_counter = 0
                for z in range(len(pileup)):
                    if pileup[z] == "-":
                        gap_counter = gap_counter + 1
                if vector[ii - 15] == ["-", "-"] and gap_counter >= len(pileup) * 0.8:
                    temp_list.append(0)
                    continue

                # NaN-Checking: case in vector part ii is no match with any element from pile-up-list
                nan_controll = 0
                for sub in range(len(pileup)):
                    if vector[ii - 15][0] == pileup[sub] or vector[ii - 15][1] == pileup[sub]:
                        nan_controll = 1

                if nan_controll == 0:
                    temp_list.append("NaN")
                    continue

                # Loop over vector
                # set temp_value to zero for every ii-loop
                temp_value = 0
                for subi in range(len(pileup)):

                    if vector[ii - 15][0] == pileup[i] and vector[ii - 15][1] == pileup[i]:
                        #   log10(1 - 10 ^ (-qscore[subi] / 10))
                        temp_pow = math.pow(10, (-pileup_qual[subi] / 10))
                        temp = math.log10(1 - temp_pow)
                        temp_value = temp_value + temp

                    elif pileup[subi] != vector[ii - 15][0] and pileup[subi] != vector[ii - 15][1]:
                        #   -qscore[subi] / 10 +log10(0.25)
                        temp = (-pileup_qual[subi] / 10 + math.log10(0.25))
                        temp_value = temp_value + temp

                    else:
                        #   log10(0.5*(1-10^(-qscore[subi]/10)) + 0.125*10^(-qsccore[subi] / 10))
                        temp_pow = math.pow(10, (-pileup_qual[subi] / 10))
                        temp = 0.5 * (1 - temp_pow) + 0.125 * temp_pow
                        temp_log = math.log10(temp)
                        temp_value = temp_value + temp_log

                temp_list.append(temp_value)
            ematrix.append(temp_list)
            temp_list = []

        # case: if reference-base at i is not a gap
        else:
            # get vector:
            if upd_reference[i] == "A":
                vector = vector_A
            elif upd_reference[i] == "C":
                vector = vector_C
            elif upd_reference[i] == "G":
                vector = vector_G
            elif upd_reference[i] == "T":
                vector = vector_T
            else:
                print("Critical Error at creating emission-matrix!")

            for ii in range(15):
                temp_value = 0

                # Checking: case: if Vector at Gap, Gap and >80% of reads are gaps:
                gap_counter = 0
                for z in range(len(pileup)):
                    if pileup[z] == "-":
                        gap_counter = gap_counter + 1
                if vector[ii] == ["-", "-"] and gap_counter >= len(pileup) * 0.8:
                    temp_list.append(0)
                    continue

                # NaN-Checking: case in vector part ii is no match with any element from pile-up-list
                nan_controll = 0
                for sub in range(len(pileup)):
                    if vector[ii][0] == pileup[sub] or vector[ii][1] == pileup[sub]:
                        nan_controll = 1

                if nan_controll == 0:
                    temp_list.append("NaN")
                    continue

                # Loop over vector
                for subi in range(len(pileup)):

                    if vector[ii - 15][0] == pileup[i] and vector[ii - 15][1] == pileup[i]:
                        #   log10(1 - 10 ^ (-qscore[subi] / 10))
                        temp_pow = math.pow(10, (-pileup_qual[subi] / 10))
                        temp = math.log10(1 - temp_pow)
                        temp_value = temp_value + temp

                    elif pileup[subi] != vector[ii - 15][0] and pileup[subi] != vector[ii - 15][1]:
                        #   -qscore[subi] / 10 +log10(0.25)
                        temp = (-pileup_qual[subi] / 10 + math.log10(0.25))
                        temp_value = temp_value + temp

                    else:
                        #   log10(0.5*(1-10^(-qscore[subi]/10)) + 0.125*10^(-qsccore[subi] / 10))
                        temp_pow = math.pow(10, (-pileup_qual[subi] / 10))
                        temp = 0.5 * (1 - temp_pow) + 0.125 * temp_pow
                        temp_log = math.log10(temp)
                        temp_value = temp_value + temp_log

                temp_list.append(temp_value)

            for ii in range(15, 30):
                temp_list.append("NaN")
            ematrix.append(temp_list)
            temp_list = []

    return ematrix


# Test-Values for testing emission-matrix
updated_sam = ""
updated_reference = "A-GN-"
build_emissionmatrix(updated_sam, updated_reference)


# Test-Values for testing emission-matrix
updated_sam = ""
updated_reference = "A-GN-"
test_ematrix = build_emissionmatrix(updated_sam, updated_reference)
print("Ematrix ready.")
