# encoding: utf-8

import click
from plasmidpicker.core import MinHash
from Bio import SeqIO
from glob import glob
from tqdm import tqdm
import h5py

@click.command()
@click.option("-i", "--indir", default=".", help="Specify FASTA input file [.]")
@click.option("-k", "--kmer", default=21, type=int, help="Threshold of probability (0 < k) [7]")
@click.option("-s", "--sketch", default=1000, type=int, help="Threshold of probability (0 < k) [7]")
@click.option("-o", "--outfile", default="plasmids_hash.h5", help="Threshold of probability [plasmids_hash.h5]")
def cmd(indir, kmer, sketch, outfile):
    fna_files = glob(indir + "/*.fna")
    fasta_files = glob(indir + "/*.fasta")
    infiles = fna_files + fasta_files

    minhash = MinHash(sketch, kmer)

    with h5py.File(outfile, "w") as f:
        pbar = tqdm(infiles)
        file_format = "fasta"
        for index, infile in enumerate(pbar):
            pbar.set_description(infile)
            record = SeqIO.read(infile, file_format)
            minhash_array = minhash.get_minhash_value(str(record.seq))

            f.create_dataset(record.id, data=minhash_array)

def main():
    cmd()

if __name__ == "__main__":
    main()
