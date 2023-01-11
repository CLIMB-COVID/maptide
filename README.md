# `basemap`

## Setup
Current setup instructions require Rust to be installed.

Installation instructions for Rust can be found here: https://www.rust-lang.org/tools/install

Once Rust is installed:
```
$ git clone https://github.com/CLIMB-COVID/basemap.git
$ cd basemap/
$ python -m venv env
$ source env/bin/activate
$ pip install --upgrade pip
$ pip install .
```

## Usage
##### Frequencies over all positions in the reference:
```
$ basemap all /path/to/file.bam
```

All arguments:

```
$ basemap all -h
usage: basemap all [-h] [--mapq MAPQ] [--baseq BASEQ] [--stats] [--decimals DECIMALS] bam

positional arguments:
  bam                  Path to BAM file

optional arguments:
  -h, --help           show this help message and exit
  --mapq MAPQ          Minimum mapping quality (default: 0)
  --baseq BASEQ        Minimum base quality (default: 0)
  --stats              Output additional per-position statistics
  --decimals DECIMALS  Number of decimal places to display (default: 3)
```

##### Frequencies over a specific region (with an index file):
```
$ basemap iquery /path/to/file.bam chrom:start-end --index /path/to/anotherfile.bai
```

All arguments:

```
usage: basemap iquery [-h] [--index INDEX] [--mapq MAPQ] [--baseq BASEQ] [--stats] [--decimals DECIMALS] bam region

positional arguments:
  bam                  Path to BAM file
  region               Region to view, specified in the form CHROM:START-END

optional arguments:
  -h, --help           show this help message and exit
  --index INDEX        Path to index (BAI) file (default: </path/to/bam>.bai)
  --mapq MAPQ          Minimum mapping quality (default: 0)
  --baseq BASEQ        Minimum base quality (default: 0)
  --stats              Output additional per-position statistics
  --decimals DECIMALS  Number of decimal places to display (default: 3)
```

##### Frequencies over a specific region (without an index file):
```
$ basemap query /path/to/file.bam chrom:start-end
```

All arguments:

```
usage: basemap query [-h] [--mapq MAPQ] [--baseq BASEQ] [--stats] [--decimals DECIMALS] bam region

positional arguments:
  bam                  Path to BAM file
  region               Region to view, specified in the form CHROM:START-END

optional arguments:
  -h, --help           show this help message and exit
  --mapq MAPQ          Minimum mapping quality (default: 0)
  --baseq BASEQ        Minimum base quality (default: 0)
  --stats              Output additional per-position statistics
  --decimals DECIMALS  Number of decimal places to display (default: 3)
```
