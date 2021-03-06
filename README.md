![PlasmidPicker](./logo.png)

## Requirements
- numpy
- cython
- biopython
- click
- mmh3
- h5py

## Installation

```
git clone https://github.com/haradama/PlasmidPicker.git
cd PlasmidPicker
pip install -r requirements.txt
python setup.py install
```

## Usage

```
Usage: plasmidpicker pick [OPTIONS]

  Detect plasmid sequence data in metagenome.

Options:
  -i, --infile TEXT        Specify FASTA input file
  -t, --threshold INTEGER  Threshold of probability (0 <= t <= 100) [70]
  -l, --length INTEGER     Threshold of contig's length [1000]
  -d, --outdir TEXT        Output directory [.]
  --help                   Show this message and exit.
```

```
Usage: plasmidpicker classify [OPTIONS]

  Perform similarity search of plasmids using MinHash.

Options:
  -i, --infile TEXT     Specify FASTA input file
  -s, --sketch INTEGER  Sketch size (1 <= s <=1000) [1000]
  -l, --length INTEGER  Threshold of contig's length [1000]
  -hi, --hits INTEGER   Number of hits to display [5]
  -o, --outfile TEXT    Output file [result.txt]
  --help                Show this message and exit.
```

## License
