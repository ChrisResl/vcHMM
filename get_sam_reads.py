import pysam


def get_sam():
    """
    This code reads a sam-file and returns an dictionary with readnames as keys and a list with startposition, readsequence, cigarstring, queryqualities.
    return: All reads from given sam-file.
    """
    samfile = pysam.AlignmentFile("data/example.sam", "rb")
    samdict = {}
    for read in samfile:
        samdict[read.query_name] = [read.reference_start,
                                    read.query_sequence, read.cigarstring, read.mapping_quality, read.query_qualities]
    return samdict


def main():
    dict = get_sam()

    for key, value in dict.items():
        print(key, ':', value)


if __name__ == '__main__':
    main()
