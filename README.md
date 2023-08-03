# vcf_parse.py
Tool for parsing VCF files and making a text report.

## Install

Install Miniconda from https://conda.io/miniconda.html

A virtual environment called vcf_parse can then be created using the following command:

```

conda env create -f env/vcf_parse.yaml

source activate vcf_parse

```

## Run


```
usage: vcf_parse.py [-h] [-v] [-i INPUT VCF] [-O OUTPUT]
                    [-t TRANSCRIPTS] [-n NTC]
                    
summary:
Takes a VCF file and parses the variants to produce a tab delimited
variant report.

positional arguments:
  input                 Filepath to input VCF file. REQUIRED.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

  -i INPUT VCF --inp INPUT VCF
                        VCF file and filepath

 -o OUTPUT, --outfile OUTPUT
                        Filepath to folder where output reports will be saved.
                        If missing, defaults to current directory.

-t PREFERRED TRANSCRIPT --pref_trans PREFERRED TRANSCRIPT
                        Preferred transcript file and filepath   
```

### Testing
To run unit tests run `python -m unittest test`
