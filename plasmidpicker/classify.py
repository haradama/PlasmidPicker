# coding: utf-8

from plasmidpicker.core  import MinHash
from Bio import SeqIO
import numpy as np
import h5py

class Classify(MinHash):
    def __init__(self, num=1000, k_length=16):
        MinHash.__init__(num, k_length)
        db_path="../data/plasmidDB.h5"
        self.num = num
        self.ref_plasmids = dict() 
        with h5py.File(db_path, "r") as h5f:
            for key in h5f.keys():
                self.ref_plasmids[key] = h5f[key].value
    
    def getSim(self, array1, array2):
        return (array1 == array2).sum() / self.num
    
    def getSimwithRef(self, seq, hits_num):
        minhash_array = self.getMin(seq)
        similarity_dict = dict()
        for key, refminhash_array in self.ref_plasmids.items():
            similarity = self.getSim(minhash_array, refminhash_array)
            if similarity:
                similarity_dict[key] = similarity
        
        for index, (key, value) in enumerate(sorted(similarity_dict.items(), key=lambda x: -x[1])):
            if index >= hits_num:
                break
            yield (key, value)
    
    def getSimandWrite(self, infile, outfile, hits_num):
        with open(outfile, "w") as fw:
            file_format = "fasta"
            for record in SeqIO.parse(infile, file_format):
                fw.write("# {0}\n".format(record.id))
                for key, value in self.getSimwithRef(str(record.seq), hits_num):
                    fw.write("{0}\t{1}\n".format(key, value))


if __name__ == "__main__":
    pass