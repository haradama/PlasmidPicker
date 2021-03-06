# encoding: utf-8

from plasmidpicker.core import CGR
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm
from datetime import datetime
import numpy as np
import json
import os

class Pick(CGR):
    def __init__(self, k_length=7):
        CGR.__init__(k_length)
        param_path = os.path.dirname(__file__) + "/data/param.json"
        with open(param_path, "r") as fr:
            param = json.loads(fr.read())
        self.weights = param["weights"]
        self.bias = param["bias"]

    def logistic_regression(self, cgr_array):
        z = np.dot(cgr_array, self.weights)
        return 1.0 / (1 + np.exp(-z)) + self.bias

    def cat_with_reverse_complement(self, seq):
        return seq + Seq("N") + seq.reverse_complement()

    def pick_plasmid(self, infile, threshold=70, length=1000, outdir=None):
        assert infile.endswith(".fasta") or infile.endswith(".fna") or infile.endswith(".fa"), "Your input file must be in .fa, .fna or .fasta format."

        nowtime = datetime.today().strftime("%Y%m%d%H%M%S")
        if outdir is None:
            outdir = os.getcwd()
        os.chdir(outdir)
        infile = os.path.abspath(infile)

        with open("plasmids{0}.fna".format(nowtime), "w") as fpl, open("chromosomes{0}.fna".format(nowtime), "w") as fch, open("unassigned{0}.fna".format(nowtime), "w") as fun:
            file_format = "fasta"
            threshold *= 0.01
            plasmid_num = 0
            chromosome_num = 0
            unassigned_num = 0

            with open(infile, "r") as f:
                total = sum([1 for line in f if line.startswith(">")])

            pbar = tqdm(SeqIO.parse(infile, file_format), total=total)
            for record in pbar:
                pbar.set_description(record.id)

                if len(record.seq) >= length:
                    self.add_seq(str(self.cat_with_reverse_complement(record.seq)))
                    cgr_array = self.get_cgr_array()
                    probability = self.logistic_regression(cgr_array.flatten())
                    record.description += " p={0}".format(str(round(probability, 3)))

                    if probability >= threshold:
                        f = fpl
                        plasmid_num += 1
                    elif probability <= (1 - threshold):
                        f = fch
                        chromosome_num += 1
                    else:
                        f = fun
                        unassigned_num += 1

                    SeqIO.write(record, f, file_format)

            print("plasmid:", plasmid_num)
            print("chromosome:", chromosome_num)
            print("unassigned:", unassigned_num)


if __name__ == "__main__":
    pass
