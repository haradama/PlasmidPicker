# encodling: utf-8

import click
from Bio.SeqRecord import SeqRecord
from glob import glob
from tqdm import tqdm
from Bio import SeqIO
import numpy as np

@click.command()
@click.option("-i", "--indir", default=".", help="Input directory [.]")
@click.option("-l", "--length", default=5000, type=int)
@click.option("-n", "--num", default=10, type=int)
@click.option("-o", "--outdir", default=".", help="Output directory [.]")
def cmd(indir, length, num, outdir):
    fna_files = glob(indir + "/*.fna")
    fasta_files = glob(indir + "/*.fasta")
    infiles = fna_files + fasta_files
    file_format = "fasta"

    if outdir.endswith("/"):
        outdir = outdir[:-1]

    pbar = tqdm(infiles)
    for infile in pbar:
        pbar.set_description(infile)
        record = SeqIO.read(infile, file_format)
        digit = len(str(num))
        if len(record.seq) - length > 0:
            positions = np.random.randint(0, len(record.seq) - length, (num))
            for index, p in enumerate(positions):
                contig = record.seq[p:p+length]
                out_record = SeqRecord(contig)
                out_record.id = record.id + "_len{0}.{1}".format(length, str(index + 1).zfill(digit))
                SeqIO.write(out_record, outdir + "/{0}.fna".format(out_record.id), file_format)

def main():
    cmd()

if __name__ == "__main__":
    main()
