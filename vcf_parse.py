#!/anaconda3/envs/python2/bin/python

"""
vcf_parse.py

Takes a VCF file and parses the variants to produce a tab delimited 
variant report.

Usage:  vcf_parse.py [-h] [-O OUTPUT] [-s SETTINGS] 
                     [-t TRANSCRIPTS] [-T TRANSCRIPT_STRICTNESS] 
                     [-b BED | -B BED_FOLDER] 
                     [-l] [-v]
                     input
        vcf_parse.py -h for full description of options.

Author:     Erik Waskiewicz
Created:    31 Aug 2018
Version:    1.0.0
Updated:    23 Oct 2018
"""
__version__ = '1.0.0'
__updated__ = '23 Oct 2018'


import argparse
import logging

from scripts.vcf_report import vcf_report
from scripts.preferred_transcripts import preferred_transcripts
from scripts.bed_object import bed_object


## -- PARSE INPUT ARGUMENTS -------------------------------------------

def get_args():
    '''
    Use argparse package to take arguments from the command line and 
    store them as an argparse object. See descriptions for full detail 
    of each argument.
    '''

    # Make empty argparse object
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=
'''
summary:
  Takes a VCF file and parses the variants to produce a tab delimited 
  variant report.
'''
    )


    # Version info
    parser.add_argument(
        '-v', '--version', action='version', 
        version=
            '%(prog)s\nversion:\t{}\nlast updated:\t{}'.format(
                __version__, __updated__
            )
    )


    # Arguments (see help string for full descriptions):
    # REQUIRED: VCF file input
    parser.add_argument(
        'input', action='store', 
        help='Filepath to input vcf file. REQUIRED.'
    )


    # OPTIONAL: Output folder, defaults to current directory if empty
    parser.add_argument(
        '-O', '--output', action='store', 
        help=
'''
Filepath to folder where output reports will be saved. 
If missing, defaults to current directory.
\n'''
    )


    # OPTIONAL: List of preferred transcripts
    parser.add_argument(
        '-t', '--transcripts', action='store', 
        help=
'''
Filepath to preferred transcripts file. 

Must be a tab seperated file with preferred transcripts in the seond 
column. If missing, all entries in the preferred transcript column 
will be labelled as 'Unknown'.
\n'''
    )


    # OPTIONAL: List of preferred transcripts
    parser.add_argument(
        '-T', '--transcript_strictness', action='store', default='low', 
        help=
'''
Strictness of matching while annotating preferred transcripts.
Default setting is low.

Options: 

high - Transcripts must be an exact match. 
       e.g. NM_001007553.2 and NM_001007553.1 won't match,
            NM_001007553.1 and NM_001007553.1 will.

low  - Transcripts will match regardless of the version number. The 
       version number is after the . at the end of a transcript 
       e.g. NM_001007553.2 and NM_001007553.1 will match.
\n'''
    )


    # OPTIONAL: either a single BED file or a folder containing BED 
    # files, only one of these can be used
    bed_files = parser.add_mutually_exclusive_group()
    bed_files.add_argument(
        '-b', '--bed', action='store', 
        help=
'''
Filepath to a single BED file. 
Cannot be used together with -B flag.
\n'''
    )
    bed_files.add_argument(
        '-B', '--bed_folder', action='store', 
        help=
'''
Filepath to folder of BED files. 
Cannot be used together with -b flag.
\n'''
    )


    # OPTIONAL: File containing the headers for the report
    parser.add_argument(
        '-s', '--settings', action='store', 
        help=
'''
Filepath to settings file. 

This is a tab seperated text file containing a number of rows, where 
each row specifies an annotation to be included in the variant report.
Only annotations included in the settings file will be included in the
variant report.
The columns in the variant report will be in the same order as the order
in which the annotations appear in the settings file.

Each row contains two columns:

Column 1 - annotations headers, these must match up with how they appear
           in the VCF (case sensitive), 

Column 2 - location where to find the data within the VCF, used to select 
           the correct parsing function.
           Either info, format, vep, filter or pref.

To make a settings file with all available options from a VCF, run:
    vcf_parse -l path_to_input_vcf > settings.txt
\n'''
    )


    # OPTIONAL: Lists all headers in a vcf then exits
    parser.add_argument(
        '-l', '--settings_list', action='store_true', 
        help=
'''
Return a list of all availabile settings to the screen, then exit.
See settings flag for usage.
\n'''
    )


    return parser.parse_args()


# -- CALL FUNCTIONS ---------------------------------------------------

if __name__ == '__main__':
    # add logger
    logger = logging.getLogger('vcf_parse')
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '%(levelname)s\t%(asctime)s\t%(name)s\t%(message)s'
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    print('{}'.format('---'*30))
    logger.info('running vcf_parse.py...')


    # Load arguments,make vcf report object and load data
    args = get_args()
    report = vcf_report()
    report.load_data(args.input, args.output)


    # If -l flag called, print headers and exit
    if args.settings_list:
        report.list_settings()
        exit()


    # If settings file provided, load settings
    if args.settings:
        report.settings(args.settings)
    else:
        logger.info('no settings file found -- outputting all data from VCF.')


    # Make variant report of whole VCF
    report.make_report()


    # If preferred transcripts provided, apply preferred transcripts to
    # variant report
    if args.transcripts:
        pt = preferred_transcripts()
        pt.load(args.transcripts)
        pt.apply(report.report_path, args.transcript_strictness)
    else:
        logger.info('no preferred transcripts file provided -- preferred ' +
        'transcripts column will all be labelled as "Unknown"')


    # If single BED file provided, make variant report with BED file 
    # applied
    if args.bed:
        bed = bed_object()
        bed.apply_single(args.bed, report)

    # If folder of BED file provided, make a seperate variant report 
    # for each BED file. Output will be saved in a folder named the 
    # same as the BED file folder, within the output directory.
    elif args.bed_folder:
        bed = bed_object()
        bed.apply_multiple(args.bed_folder, report)

    else:
        logger.info('no BED files provided')

    # Finish
    logger.info('vcf_parse.py completed\n{}'.format('---'*30))
