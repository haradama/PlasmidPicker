# encoding: utf-8

import click
from plasmidpicker.core import CGR
from Bio import SeqIO
import numpy as np
from glob import glob
from tqdm import tqdm

@click.command()
@click.option("-i", "--indir", default=".", help="Intput directory [.]")
@click.option("-k", "--kmer", default=7, type=int)
@click.option("-o", "--outfile", default=".", help="Output directory [.]")
def cmd(indir, kmer, outfile):
    fna_files = glob(indir + "/*.fna")
    fasta_files = glob(indir + "/*.fasta")

    infiles = fna_files + fasta_files

    cgr = CGR(kmer)
    pbar = tqdm(infiles)
    file_format = "fasta"
    array_size = int((4**kmer)**0.5)
    for index, infile in enumerate(pbar):
        pbar.set_description(infile)
        record = SeqIO.read(infile, file_format)

        cgr.add_seq(str(seq + Seq("N") + seq.reverse_complement()))
        cgr_array = cgr.get_cgr_array()

        if index == 0:
            cgr_arrays = cgr_array.reshape(1, array_size, array_size)

        else:
            cgr_array = cgr_array.reshape(1, array_size, array_size)
            cgr_arrays = np.append(cgr_arrays, cgr_array, axis=0)

    np.save(outfile, cgr_arrays)

def main():
    cmd()

if __name__ == "__main__":
    main()
