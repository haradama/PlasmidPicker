# coding: utf-8

from keras import models
from Bio import SeqUtils
from Bio import SeqIO
import pandas as pd
from cgr import CGR
import numpy as np
import sys
import os

args = sys.argv
input_fasta = args[1]
output_dir = args[2]

if output_dir.endswith("/"):
    output_dir = output_dir.rstrip("/")

if output_dir not in os.listdir("."):
    os.mkdir(output_dir)


with open("./harada_MLP2.json") as f:
    json = f.read()

cgr_array = []
gc_content_array = []
for record in SeqIO.parse(input_fasta, "fasta"):
    make_cgr = CGR(record.seq, 7)
    cgr = make_cgr.make_representation()
    cgr_array.append(list(cgr))
    gc_content = SeqUtils.GC(record.seq)
    gc_content_array.append(gc_content * 0.01)

model = models.model_from_json(json)
model.load_weights('harada_MLP2.hdf5')

print("model is loaded.")

cgr_array = np.array(cgr_array)
x, y, _ = cgr_array.shape
cgr_array = cgr_array.reshape((x, 1, y, y))
gc_content_array = np.array(gc_content_array)

input_data = [cgr_array, gc_content_array]

results = model.predict(input_data)

label = ["chromosome", "plasmid"]
plasmids = []
chromosome = []
matrix = []
for result, record in zip(results, SeqIO.parse(input_fasta, "fasta")):
    argmax = np.argmax(result)
    row = [record.id]
    row.extend(result)
    if argmax:
        plasmids.append(record)
        row.append(label[argmax])
    else:
        chromosome.append(record)
        row.append(label[argmax])

    matrix.append(row)

df = pd.DataFrame(matrix, columns=["id", "prob_chromosme", "prob_plasmid", "predict"])
df.to_csv(output_dir + "/" + "log.csv")
SeqIO.write(plasmids, output_dir + "/" + "plasmids.fna", "fasta")
SeqIO.write(chromosome, output_dir + "/" + "chromosome.fna", "fasta")

print("finished")