"""
This program allows to align reads with a reference sequence by the Burrows-Wheeler transform.

Authors : Pichon Julien , Rambaud Opale

"""





from operator import itemgetter
from collections import Counter
import copy
import sys
import argparse
import os
import csv


    
def isfile(path):
    """
    Check the existence of the file 
    """

    if not os.path.isfile(path):

        if os.path.isdir(path):
            err = f"{path} is a directory"
        else:
            err = f"{path} does not exist"

        raise argparse.ArgumentTypeError(err)

    return path


def arguments():
    """
    arguments setting
    """

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-f", "--ref_file", dest = "ref_file",
        type = isfile, required = True)

    parser.add_argument("-r", "--reads_file", dest = "reads_file",
        type = isfile, required = True,)

    args = parser.parse_args()

    return args.ref_file, args.reads_file

def read_ref_sequence(fasta_file):
    """
    Return a reference sequence from initial fasta file
    """

    with open(fasta_file, "r") as filin:
        ref_sequence = ""
        for line in filin:
            if not line.startswith(">"):
                ref_sequence += line.strip()

    return ref_sequence


def read_reads(fasta_file):
    """
    Return the reads from a reads file in fasta format 
    """

    reads_dict = {}

    with open(fasta_file, "r") as f_in:
        for line in f_in:
            if line.startswith(">"):
                read_id = line[1:].strip()
                read_id = read_id.split(".")[0].split("-")[0].split("+")[0]
                read_id = read_id[:-1]
            else:
                reads_dict[read_id] = line.strip()

    return reads_dict


def bwt(seq):   
    """
    Return a sequence transformed by the Burrows-Wheeler transform 
    """

    seq = "$" + seq
    seq_l = len(seq)
    
    seq_sorted = sorted([seq[i:seq_l] + 
        seq[0:i] for i in range(1, seq_l + 1)])
    
    index_seq = seq_sorted.index(seq)
    transformed_seq = ''.join([j[-1] for j in seq_sorted])

    return index_seq, transformed_seq


def bwt_restore(index_seq, transformed_seq):
    """
    Restore the reference sequence
    """

    transformed_l, index_transformed = get_position(index_seq, transformed_seq)
    index_seq_init = [index_seq]

    for i in range(1, transformed_l - 1):
        index_seq_init.append(index_transformed[index_seq_init[i - 1]])
    seq_init_reversed = [transformed_seq[i] for i in index_seq_init]
    seq_init_reversed.reverse()
    restored_seq = ''.join(seq_init_reversed)
    
    return restored_seq


def build_tally(transformed_seq):
    """
    Generate a tally table based on the seuence transformed by the Burrows-Wheeler transform
    """

    dico, tally_table = {"$": 0, "A" : 0, "C" : 0, "G" : 0, "T" : 0}, []
    
    for n in transformed_seq:
        dico[n] += 1
        tally_table.append(copy.deepcopy(dico))

    return tally_table


def get_count(transformed_seq):
    """
    Count the characters and return a dictionnary 
    """

    nucl, count = {'$': 0, 'A': 0, 'C': 0, 'G': 0, 'T': 0}, Counter(sorted(transformed_seq))
    nb = 0
    for n in nucl.keys():
        nucl[n] = nb
        nb += count[n]

    return nucl


def get_position(index_seq, transformed_seq):
    """
    Return pos of each character in the initial sequence
    """

    transformed_l = len(transformed_seq)
    
    couple_index_base = sorted([(i, x) for i, x in enumerate(transformed_seq)], key = itemgetter(1))
    index_transformed = [None for i in range(transformed_l)]
    
    for index, couple in enumerate(couple_index_base):
        j,_ = couple
        index_transformed[j] = index

    return transformed_l, index_transformed


def mapping(read, infs, tally_table):
    """
    Map the read
    """

    a = len(read) - 1
    inf = infs[read[a]]
    supp = infs[read[a]] + tally_table[-1][read[a]] - 1

    while a >= 1 and inf <= supp:
        a -= 1
        inf = infs[read[a]] + tally_table[inf - 1][read[a]]
        supp = infs[read[a]] + tally_table[supp][read[a]] - 1

    return inf, supp


def display_alignment(ref, read, index_transformed, marke):
    """
    Return the alignment of a read with the reference sequence
    """

    poss = []
    for i in range(marke[0], marke[1] + 1):
        poss.append(index_transformed[i] + 1)
    
    poss.sort()
    nb_match = len(poss)

    align = "." * len(seq_ref)
    for pos in poss:
        alig = align[:pos] + read + align[pos + len(read):]
    

    return align, nb_match


def get_matches(index_transformed, marke):
    """
    Return the number of matches for each read
    """

    pos = []
    for i in range(marke[0], marke[1] + 1):
        pos.append(index_transformed[i] + 1)
    pos.sort()
    matches = len(pos)

    return pos, matches

def fill(text, width = 80):
    """
    Format the text in the good fasta format 
    """

    return os.linesep.join(text[i:i + width] for i in range(0, len(text), width))



def write_transform(transformed_seq):
    """
    Whrite the result of the Burrows-Wheeler transform in a text file
    """

    with open("transformed_seq.txt", "w") as filout:
        filout.write("> Burrows-Wheeler Transform" + "\n")
        filout.write(fill(transformed_seq))


def write_align(align):
    """
    Generate a CSV file containing : read's ID, number of matches
    and pos for each match
    """

    with open("align.csv", "w") as filout:
        fields = ["Read", "Matches", "pos"]
        f_writer = csv.DictWriter(filout, fieldnames = fields)
        f_writer.writeheader()
        for read, value in align.items():
            content = {"Read" : read,
                    "Matches": value["Matches"],
                    "pos" : value["pos"]}
            f_writer.writerow(content)
            
def main():
    """
    Main program function
    """

    #1 : set the arguments
    ref_file, reads_file = arguments()

    #2 : obtain the ref seq and the read parsing 
    ref_sequence = read_ref_sequence(ref_file)
    reads_dict = read_reads(reads_file)
    
    #3 : apply the BWT on the ref seq
    index_seq, transformed_seq = bwt(ref_sequence)
    
    #4 : restore the ref seq 
    rest_seq = bwt_restore(index_seq, transformed_seq)
    
    
    #5: Build the tally table
    tally_table = build_tally(transformed_seq)

    #6 : count character occurence
    infs = get_count(transformed_seq)

    #7 : reference position of each characters 
    transformed_l, index_transformed = get_position(index_seq, transformed_seq)

    #8 : get the alignments
    align = {}  
    for id_read, read in reads_dict.items():
        
        marke = mapping(read, infs, tally_table)
        pos, matches = get_matches(index_transformed, marke)
        align[id_read] = {"Matches" : matches, "pos" : pos}

    #9 : write results in a csv file for alignment and in a text file for transformd 
    write_align(align)
    write_transform(transformed_seq)
            


if __name__ == "__main__":
    main()
