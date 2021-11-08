"""
This file was taken from https://github.com/ksahlin/isONcorrect/blob/master/scripts/exon_experiment/simulate_reads.py and has been altered after consulting the author
Run with: python generateTestCasesSimple.py  --outfolder testout


"""

import os, sys
import random
import argparse
import errno
import math


'''
    Below code taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''

def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
#converts the list of characters (used for the dynamic sequence) into a string (immutable)
def convert(s):
    # initialization of string to ""
    new = ""

    # traverse in the string
    for x in s:
        new += x

        # return string
    return new
#generates a random sequence with length 'length'
def generate_random_sequence_by_length(length):
    sequence=[]
    for j in range(0, length):
        base = random.randrange(3)
        if base == 0:
            sequence.append("A")
        elif base == 1:
            sequence.append("C")
        elif base == 2:
            sequence.append("G")
        elif base == 3:
            sequence.append("T")
    final_sequence = convert(sequence)
    return final_sequence


#checks whether the given set of arguments matches the expected ones
def check_valid_args(args, ref):
    # assert args.start < args.stop and args.start < min(args.coords)
    # assert args.stop > args.start and args.stop > max(args.coords)
    assert args.coords == sorted(args.coords)  # has to be in increasing order
    assert len(ref[list(ref.keys())[0]][0]) >= max(args.coords)
    print(args.coords, args.probs)
    assert (len(args.coords) - 1) == len(args.probs)


#used to add errors to a read
def simulate_reads(isoforms, error_lvls):
    reads = {}

    for i_acc, isoform in isoforms.items():
        read = []
        qual = []
        was_del = False
        #iterate through the sequence
        if not i_acc=="@sim|correct|full":
            for l, n in enumerate(isoform):

                p_correct_reading = error_lvls
                p_error = 1.0 - p_correct_reading
                r = random.uniform(0, 1)
                if r > p_correct_reading:
                    error = True
                else:
                    error = False

                if error:
                    r = random.uniform(0, 1)
                    if r < 0.4:  # deletion(those values depend on current base caller )
                        was_del = p_error
                        pass
                    elif 0.4 <= r < 0.7: #substitution
                    #we do not want to substitute the same base, so we drop the current base from sub_possibilities
                        sub_possibilities="ACGT".replace(n,'')
                        read.append(random.choice(sub_possibilities))
                        qual.append(round(-math.log(p_error, 10) * 10))

                    else: #insertion
                        read.append(n)
                        qual.append(round(-math.log(p_error, 10) * 10))

                        r_ins = random.uniform(0, 1)
                        ins_len=1
                        while r_ins >= 0.7:
                            ins_len += 1
                            read.append(random.choice("ACGT"))
                            r_ins = random.uniform(0, 1)
                            qual.append(round(-math.log(0.7, 10) * 10))

                else:
                    if was_del:  # adding uncertainty from previous deleted base
                        read.append(n)
                        qual.append(round(-math.log(was_del, 10) * 10))
                    else:
                        read.append(n)
                        qual.append(round(-math.log(p_error, 10) * 10))
                    was_del = False
            if not read:
                continue
            read_seq = "".join([n for n in read])
            qual_seq = "".join([chr(q + 33) for q in qual])
            reads[i_acc]=(read_seq, qual_seq)
        else:
            read_seq=isoform
            qual_seq = "".join([chr(q + 33) for q in qual])
            reads[i_acc] = (read_seq, qual_seq)

    return(reads)