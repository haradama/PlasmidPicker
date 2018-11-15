# coding: utf-8

import click
from plasmidpicker.pick import Pick
from plasmidpicker.classify import Classify

@click.group()
def cmd():
    pass

@cmd.command(help="Detect plasmid sequence data in metagenome.")
@click.option("-i", "--infile", help="Specify FASTA input file")
@click.option("-t", "--threshold", default=70, type=int, help="Threshold of probability (0 <= t <= 100) [70]")
@click.option("-l", "--length", default=1000, type=int, help="Threshold of contig's length [1000]")
@click.option("-d", "--outdir", help="Output directory [.]")
def pick(infile, threshold, length, outdir=None):
    if infile is None:
        raise click.BadParameter("Please specify FASTA input file")
    if threshold < 0 or threshold > 100:
        raise click.BadParameter("Threshold of probability (0 <= t <= 100) [70]")
    picker = Pick()
    picker.pick_plasmid((infile, threshold, length, outdir)

@cmd.command(help="Perform similarity search of plasmids using MinHash.")
@click.option("-i", "--infile", help="Specify FASTA input file")
@click.option("-s", "--sketch", default=1000, type=int, help="Sketch size (1 <= s <=1000) [1000]")
@click.option("-l", "--length", default=1000, type=int, help="Threshold of contig's length [1000]")
@click.option("-hi", "--hits", default=5, type=int, help="Number of hits to display [5]")
@click.option("-o", "--outfile", default="result.txt", help="Output file [result.txt]")
def classify(infile, sketch, outfile, length, hits):
    if infile is None:
        raise click.BadParameter("Please specify FASTA input file")
    if sketch < 1 or sketch > 1000:
        raise click.BadParameter("Sketch size (1 <= s <=1000) [1000]")
    classifier = Classify(sketch)
    classifier.output_similarity(infile, outfile, hits)

def main():
    cmd()

if __name__ == "__main__":
    main()
