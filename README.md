# `maptide`

## Setup
Current setup instructions require Rust to be installed.

Installation instructions for Rust can be found here: https://www.rust-lang.org/tools/install

Once Rust is installed:
```
$ git clone https://github.com/CLIMB-COVID/maptide.git
$ cd maptide/
$ python -m venv env
$ source env/bin/activate
$ pip install --upgrade pip
$ pip install .
```

## Usage
```
$ maptide -h
usage: maptide [-h] [--region REGION] [--index INDEX] [--mapq MAPQ] [--baseq BASEQ] [--noindex] [--stats] [--decimals DECIMALS] bam

positional arguments:
  bam                  Path to BAM file

optional arguments:
  -h, --help           show this help message and exit
  --region REGION      Region to view, specified in the form CHROM:START-END (default: everything)
  --index INDEX        Path to index (BAI) file (default: </path/to/bam>.bai)
  --mapq MAPQ          Minimum mapping quality (default: 0)
  --baseq BASEQ        Minimum base quality (default: 0)
  --noindex            Do not use an index file when querying the BAM file (default: False)
  --stats              Output additional per-position statistics (default: False)
  --decimals DECIMALS  Number of decimal places to display (default: 3)
```

##### Frequencies over all positions in the reference:
```
$ maptide /path/to/file.bam
```

##### Frequencies over a specific region (with an index file):
If the index file has the same path as the BAM file, but with `.bai` appended on the end: 
```
$ maptide /path/to/file.bam --region chrom:start-end
```

Otherwise, the path needs to be specified:
```
$ maptide /path/to/file.bam --region chrom:start-end --index /path/to/index.bai
```

##### Frequencies over a specific region (without an index file):
```
$ maptide /path/to/file.bam chrom:start-end --region chrom:start-end --noindex
```
