def find_base_state(xtilde, upd_ref):
    """
    This code finds the base state (Hidden States) at every position of updated reference.
    e.g. at R(i = 7) base is gap, and there are only Gs within the reads: -> Base State is [G, G]
    :param xtilde:  Hidden State at position Ri
    :param upd_ref: Reference Genome with gaps
    :return:
    """
    base_state = []

    #Builds Vectors:
    vector_A = [["A", "A"], ["C", "C"], ["G", "G"], ["T", "T"], ["A", "C"], ["A", "G"], ["A", "T"], ["A", "-"],
                ["C", "G"], ["C", "T"], ["C", "-"], ["G", "T"], ["G", "-"], ["T", "-"], ["-", "-"]]
    vector_C = [["C", "C"], ["A", "A"], ["G", "G"], ["T", "T"], ["A", "C"], ["C", "G"], ["C", "T"], ["C", "-"],
                ["A", "G"], ["A", "T"], ["A", "-"], ["G", "T"], ["G", "-"], ["T", "-"], ["-", "-"]]
    vector_G = [["G", "G"], ["A", "A"], ["C", "C"], ["T", "T"], ["A", "G"], ["C", "G"], ["G", "T"], ["G", "-"],
                ["A", "C"], ["A", "T"], ["A", "-"], ["C", "T"], ["C", "-"], ["T", "-"], ["-", "-"]]
    vector_T = [["T", "T"], ["A", "A"], ["C", "C"], ["G", "G"], ["A", "T"], ["C", "T"], ["G", "T"], ["T", "-"],
                ["A", "C"], ["A", "G"], ["A", "-"], ["C", "G"], ["C", "-"], ["G", "-"], ["-", "-"]]

    # Loop over updated reference:
    for i in range(len(upd_ref)):
        # Case: xtilde(i) == -1: skip this, just base = base
        if xtilde[i] == -1:
            base_state.append(-1)

        else:
            # Get Base at Ri:
            if upd_ref[i] == "A":
                vector = vector_A
                vector.extend(vector_A)
                base_state.append(vector[xtilde[i]])

            elif upd_ref[i] == "C":
                vector = vector_C
                vector.extend(vector_A)
                base_state.append(vector[xtilde[i]])

            elif upd_ref[i] == "G":
                vector = vector_G
                vector.extend(vector_A)
                base_state.append(vector[xtilde[i]])

            elif upd_ref[i] == "T":
                vector = vector_T
                vector.extend(vector_A)
                base_state.append(vector[xtilde[i]])

            elif upd_ref[i] == "N":
                base_state.append(upd_ref[i])

            elif upd_ref[i] == "-":
                vector = vector_A
                vector.extend(vector_A)
                base_state.append(vector[xtilde[i]])

            else:
                print("Critical Error at def: find_base_state. Unknown Base at Position ", i, " in updated Reference!")
                base_state.append(upd_ref[i])

    return base_state
