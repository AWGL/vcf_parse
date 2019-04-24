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
usage: vcf_parse.py [-h] [-v] [-O OUTPUT] [-t TRANSCRIPTS]
                    [-T TRANSCRIPT_STRICTNESS] [-b BED | -B BED_FOLDER]
                    [-k KNOWN_VARIANTS] [-c CONFIG] [-l] [-F]
                    input

summary:
Takes a VCF file and parses the variants to produce a tab delimited
variant report.

positional arguments:
  input                 Filepath to input VCF file. REQUIRED.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -O OUTPUT, --output OUTPUT

                        Filepath to folder where output reports will be saved.
                        If missing, defaults to current directory.

  -t TRANSCRIPTS, --transcripts TRANSCRIPTS

                        Filepath to preferred transcripts file.

                        Must be a tab seperated file with preferred transcripts in the second
                        column. If missing, all entries in the preferred transcript column
                        will be labelled as 'Unknown'.

  -T TRANSCRIPT_STRICTNESS, --transcript_strictness TRANSCRIPT_STRICTNESS

                        Strictness of matching while annotating preferred transcripts.
                        Default setting is low.

                        Options:

                        high - Transcripts must be an exact match.
                               e.g. NM_001007553.2 and NM_001007553.1 won't match,
                                    NM_001007553.1 and NM_001007553.1 will.

                        low  - Transcripts will match regardless of the version number. The
                               version number is after the . at the end of a transcript
                               e.g. NM_001007553.2 and NM_001007553.1 will match.

  -b BED, --bed BED
                        Filepath to a single BED file.

                        The BED file will be applied to the variant report and a seperate
                        report saved with the BED file applied. This report will be saved in
                        the same output folder as the original variant report, with the BED
                        file name added to it.
                        Cannot be used together with -B flag.

  -B BED_FOLDER, --bed_folder BED_FOLDER

                        Filepath to folder containing BED files.

                        Each BED file will be applied to the variant report and a seperate
                        report saved with the BED file applied. These reports will be saved in
                        a new folder within the output folder, named the same as the input BED
                        folder.
                        The file names will be the same as the original variant report, with
                        the BED file name added to them.
                        Cannot be used together with -b flag.

  -k KNOWN_VARIANTS, --known_variants KNOWN_VARIANTS

                        Filepath to known variants file.

                        This is a VCF file containing any known variants and an associated
                        classification. The classification will be added to the variant
                        report. The VCF must have an annotation named 'Classification' within
                        the INFO field for each variant.

                        Key:
                        0 - Artifact
                        1 - Benign
                        2 - Likely benign
                        3 - VUS
                        4 - Likely pathogenic
                        5 - Pathogenic

  -c CONFIG, --config CONFIG

                        Filepath to config file.

                        This is a tab seperated text file containing a number of rows, where
                        each row specifies an annotation to be included in the variant report.
                        Only annotations included in the config file will be included in the
                        variant report.
                        The columns in the variant report will be in the same order as the
                        order in which the annotations appear in the config file.

                        Each row contains:

                        Column 1 - Required. Annotation headers, these must match up with how
                                   they appear in the VCF (case sensitive).

                        Column 2 - Required. Location where to find the data within the VCF,
                                   used to select the correct parsing function.
                                   options: info, format, vep, filter or pref.

                        Column 3 - Optional. Alternative name for column header.

                        To make a config file with all available options from a VCF, run:
                            vcf_parse -l path_to_input_vcf > config.txt

  -l, --config_list
                        Return a list of all availabile config to the screen, then exit.
                        See CONFIG section for usage.

  -F, --filter_non_pass

                        Filters out any variants where the FILTER annotation is not
                        PASS. If missing then there will be no fitering based on the
                        FILTER annotation.

```
## Filtering of output

By default, from v0.1.1, there is no filtering of variants based on the filter column in the VCF.

**Important for SomaticAmplicon pipeline:**

The SomaticAmplicon pipeline requires vcf_parse to filter out any variant which had any annotation other than PASS in the filter column. This is the original behaviour in v0.1.0 and the previous vcf parsing tool.

To preserve this behavior in v0.1.1 onwards, include the `--filter_non_pass` flag.

## Special formatting for report fields:  
- % Allele frequency has been calculated from the AD
- Genotypes are reformatted from 0/0, 0/1 and 1/1 to HOM_REF, HET and HOM_ALT, respectively, to prevent them from appearing as dates in Excel
- dbSNP, COSMIC and HGMD are parsed from the 'Existing_variation' VEP field
- All population frequencies from ExAC and 1KG are re-formatted to a percentage
- Intron and exon numberings are renamed from x/y to x|y (where x and y are numbers), to prevent them appearing as dates in Excel
- HGVS coding (HGVSc) and protein (HGVSp) sequences have had the transcript name trimmed off

## Testing
To run unit tests run `python -m unittest test`
