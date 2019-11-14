"""
David McGinn
CMPS-4553
11-20-2019
This program compares a sequence indicated by the user to every other sequence
in a given folder. An example on how to do this from the command line:
python3 diff.py ./common_cold/rhinovirus_14.gb ./common_cold
samples are compared and given an alignment score. the higher the score the more
nucleotide matches and thus the more similar
"""

from Bio import SeqIO
from Bio import pairwise2
import numpy
import sys, os
from difflib import SequenceMatcher

from Bio.pairwise2 import format_alignment

def compare_samples(samp1,samp2):
    alignments = pairwise2.align.globalxx(samp1[0],samp2[0],score_only=True)
    print(samp1[1], "and",samp2[1]+"Have:",alignments,"matching nucleotides on best alignment")

if __name__ == "__main__":
    main_Seq_Path = sys.argv[1]
    main_Seq_Name = os.path.basename(main_Seq_Path).split(".")[0]

    for seq_Record in SeqIO.parse(main_Seq_Path,"genbank"):
        main_Seq = seq_Record.seq

    print(main_Seq_Name, "is", str(len(main_Seq)),"Nucleotides Long")
    print("____________________________________________________")

    folder = sys.argv[2]
    filenames = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder,f))]

    for f in filenames:
        cand_Sample_Path = folder + "/" + str(f)
        fname = f.split(".")[0]
        if fname != main_Seq_Name:
            for seq_Record in SeqIO.parse(cand_Sample_Path,"genbank"):
                cand_Sample = seq_Record.seq
            compare_samples((main_Seq,main_Seq_Name),(cand_Sample,fname))
    