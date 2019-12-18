import pysam
from Bio import SeqIO


def get_ref_fasta(file_name):
    """
    This code reads a fasta-file and return the first sequence.
    Needs: a full file-name, e.g. : ref.fa
    :return: DNA-Sequence (as String) of a given fasta-file
    """
    record_dict = list(SeqIO.parse(file_name, "fasta"))
    return record_dict[0].seq


def get_sam(samfile):
    """
    This code reads a sam-file and returns a list with startposition, readnames,
    readsequence, [cigarstring as tuples], [queryqualities].
    return: All sorted reads from given sam-file.
    """
    samfile = pysam.AlignmentFile(samfile, "rb")
    sam = []
    for read in samfile:
        sam.append([read.reference_start, read.query_name, read.query_sequence,
                    read.cigartuples, read.mapping_quality, read.query_qualities.tolist()])
        sam.sort()
    return sam


def update_reads(sam):
    """
    This code gets the cut out the hard and soft clipped entries from read sequence and/or cigar string
    and gets the insertion- and deletion positions to update reads and query qualities
    """
    insertions = []
    for read in sam:
        soft_at_beginning = True
        hard_at_beginning = True
        pos = 0
        # read[3] is cigarstring (operation, length)
        for operation in read[3]:
            # soft and hard clipping only at beginning and end
            if operation[0] == 4:
                # Soft clipped, delete sequence from read and from cigar string
                if soft_at_beginning:
                    soft_at_beginning = False
                    read = [read[0] + operation[1], read[1],
                            read[2][operation[1]:], read[3][1:], read[4], read[5]]
                else:
                    read = [read[0], read[1],
                            read[2][:-(operation[1])], read[3][:-1], read[4], read[5]]
            elif operation[0] == 5:
                # Hard clipped, delete tuple from cigar string
                if hard_at_beginning:
                    hard_at_beginning = False
                    read = [read[0], read[1], read[2],
                            read[3][1:], read[4], read[5]]
                else:
                    read = [read[0], read[1], read[2],
                            read[3][:-1], read[4], read[5]]
            elif operation[0] == 1:
                # Insertion: save position, length and readname into insertion list
                insertions.append((read[0] + pos, operation[1], read[1]))
            elif operation[0] == 2:
                # Deletion: gaps in reads einfügen und in liste speichern für folgende reads
                # insert insertion of one read into the other ones
                updated_read = read[2][: pos] + \
                    operation[1] * '-' + read[2][pos:]
                read = [read[0], read[1], updated_read,
                        read[3], read[4], read[5]]
            pos += operation[1]
        for insert in insertions:
            # insert insertiongaps into reads and into query quality
            if (read[1] != insert[2]) and (insert[0] > read[0]) and (insert[0] - read[0] < len(read[2])):
                # print('here')
                gapped_read = read[2][:insert[0] - read[0]] + \
                    insert[1] * '-' + read[2][insert[0] - read[0]:]
                # print(read[2])
                changed_qual = read[5][:insert[0] - read[0]] + \
                    insert[1] * [255] + read[5][insert[0] - read[0]:]
                # print(changed_qual)
                read = [read[0], read[1], gapped_read,
                        read[3], read[4], changed_qual]
    return sam, insertions


def del_duplicate_ins(insertions):
    """
    This function deletes all duplicate insertions from insertion list
    """
    insertions = sorted(insertions)
    unique_inserts = []
    for i in range(len(insertions) - 1):
        if insertions[i][0] != insertions[i + 1][0]:
            unique_inserts.append(insertions[i])
    return unique_inserts


def update_ref(ref_seq, insertions):
    """
    insertion of gaps into reference sequence
    """
    for ins in insertions:
        ref_seq = ref_seq[:ins[0]] + ins[1] * '-' + ref_seq[ins[0]:]
    return ref_seq


def main():
    sam = get_sam('data/example.sam')
    ref_seq = get_ref_fasta('data/ref.fa')
    updated_sam, insertions = update_reads(sam)
    unique_inserts = del_duplicate_ins(insertions)
    updated_refseq = update_ref(ref_seq, unique_inserts)
    # print(updated_refseq)
    # for read in updated_sam:
    #   print(read[5])


if __name__ == '__main__':
    main()
