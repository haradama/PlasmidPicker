# coding: utf-8

import click
from plasmidpicker.pick import Pick
from plasmidpicker.classify import Classify

@click.group()
def cmd():
    pass

@cmd.command()
@click.option('--infile', help='Prefix of greetings.')
@click.option('--threshold', default=70, type=int, help='Prefix of greetings.')
@click.option('--length', default=1000, type=int, help='Prefix of greetings.')
def pick():
    click.echo('Hello, World!')

@cmd.command()
@click.option('--infile', help='Prefix of greetings.')
@click.option('--prefix', default=1000, help='Prefix of greetings.')
@click.option('--num', default=1000, type=int, help='Prefix of greetings.')
@click.option('--output', default=70, help='Prefix of greetings.')
def classify():
    click.echo('Konnichiwa, Sekai!')

def main():
    cmd()

if __name__ == "__main__":
    main()