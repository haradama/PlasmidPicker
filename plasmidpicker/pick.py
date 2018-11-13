# encoding: utf-8

from plasmidpicker.core import CGR
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
import numpy as np
import json

class Pick(CGR):
    def __init__(self, k_length, pram_path="../data/param.json"):
        CGR.__init__(k_length)
        with open(param_path, "r") as fr:
            param = json.loads(fr.read())
        self.weights = param["weights"]
        self.bias = param["bias"]
        
    def logisticRegression(self, cgr_array):
        z = np.dot(cgr_array, self.weights)
        return 1.0 / (1 + np.exp(-z)) + self.bias
    
    def catSeqwithRCSeq(self, seq):
        return seq + Seq("N") + seq.reverse_complement()
    
    def extractPlasmidSeq(self, infile, k_length=7, threshold=0.7, seq_length_threshold=1000):
        assert infile.endswith(".fasta") or infile.endswith(".fna") or infile.endswith(".fa"), "Your input file must be in .fa, .fna or .fasta format."
        
        nowtime = datetime.today().strftime("%Y%m%d%H%M%S")
        with open("plasmidpicker_plasmids{0}.fna".format(nowtime), "w") as fpl, open("plasmidpicker_chromosomes{0}.fna".format(nowtime), "w") as fch, open("plasmidpicker_unassigned{0}.fna".format(nowtime), "w") as fun:
            file_format = "fasta"
            for record in SeqIO.parse(infile, file_format):
                if len(record.seq) >= seq_length_threshold:
                    self.add_seq(str(self.catSeqwithRCSeq(record.seq)))
                    cgr_array = self.get_cgr_array()
                    probability = self.logisticRegression(cgr_array.flatten())
                    record.description += " p={0}".format(str(round(probability, 3)))

                    if probability > threshold:
                        SeqIO.write(record, fpl, file_format)
                    elif probability < (1 - threshold):
                        SeqIO.write(record, fch, file_format)
                    else:
                        SeqIO.write(record, fun, file_format)


if __name__ == "__main__":
    pass