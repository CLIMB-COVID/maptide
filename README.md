# `maptide`

## Setup
#### Install via pip
```
$ pip install maptide
```
Depending on your operating system, the Rust compiler may need to be installed.

Installation instructions for the Rust compiler can be found here: https://www.rust-lang.org/tools/install

#### Build from source
Building from source requires the Rust compiler.

Once the Rust compiler is installed:
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
usage: maptide [-h] [-v] [-r REGION] [-i INDEX] [-m MAPPING_QUALITY] [-b BASE_QUALITY] [-s] [-d DECIMALS] bam

positional arguments:
  bam                   Path to BAM file

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -r REGION, --region REGION
                        Region to view, specified in the form CHROM:START-END (default: everything)
  -i INDEX, --index INDEX
                        Path to index (BAI) file (default: </path/to/bam>.bai)
  -m MAPPING_QUALITY, --mapping-quality MAPPING_QUALITY
                        Minimum mapping quality (default: 0)
  -b BASE_QUALITY, --base-quality BASE_QUALITY
                        Minimum base quality (default: 0)
  -s, --stats           Output additional per-position statistics (default: False)
  -d DECIMALS, --decimals DECIMALS
                        Number of decimal places to display (default: 3)
```

#### Frequencies over all positions:
```
$ maptide /path/to/file.bam
```

#### Frequencies over a region:
```
$ maptide /path/to/file.bam --region chrom:start-end
```
If a region is specified, `maptide` will check for an index file with the same path as the BAM file, but with `.bai` appended on the end (i.e. `/path/to/file.bam.bai`).

If an index file does not exist in this location, `maptide` will still run anyway, just without an index file.

Index files that do not follow the naming convention `/path/to/file.bam.bai` can still be used, but a path to the file needs to be provided:
```
$ maptide /path/to/file.bam --region chrom:start-end --index /path/to/index.bai
```
