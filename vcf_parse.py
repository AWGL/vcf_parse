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
        description='''Takes a VCF file and converts it into a tsv format 
        variant report. Can also take a BED file or a folder of BED 
        files and output a seperate variant report containing only the 
        calls within each of the BED files.'''
    )

    # Arguments (see help string for full descriptions):
    # REQUIRED: VCF file input
    parser.add_argument(
        'input', action='store', 
        help='Filepath to input vcf file. Required.'
    )

    # OPTIONAL: Output folder, defaults to current directory if empty
    parser.add_argument(
        '-O', '--output', action='store', 
        help='''Filepath to folder where output reports will be saved. If 
        empty, defaults to current directory.'''
    )

    # OPTIONAL: File containing the headers for the report
    parser.add_argument(
        '-s', '--settings', action='store', 
        help='''Filepath to settings file. This is a tab seperated file which 
        specifies the order of the report headers in the first column and the 
        location to find the data in the second column. Column headers must 
        match up with how they appear in the VCF, for a full list of the 
        available settings, run vcf_parse with the -l flag.'''
    )

    # OPTIONAL: List of preferred transcripts
    parser.add_argument(
        '-t', '--transcripts', action='store', 
        help='''Filepath to preferred transcripts file. must be a tab 
        seperated file with preferred transcripts in the seond column. If 
        empty, all entries in the preferred transcript column (if present) 
        will be false.'''
    )


    # OPTIONAL: List of preferred transcripts
    parser.add_argument(
        '-T', '--transcript_strictness', action='store', default=2, 
        help='''Strictness of matching while annotating preferred transcripts.
        Levels - 1: Transcripts must be an exact match. 2: transcripts will 
        match regardless of the version number after the . at the end of a 
        transcript (i.e. NM_001007553.2 will match with NM_001007553.1)'''
    )
    
    # OPTIONAL: either a single BED file or a folder containing BED 
    # files, only one of these can be used
    bed_files = parser.add_mutually_exclusive_group()
    bed_files.add_argument(
        '-b', '--bed', action='store', 
        help='''Filepath to a single BED file. 
        Cannot be used together with -B flag.'''
    )
    bed_files.add_argument(
        '-B', '--bed_folder', action='store', 
        help='''Filepath to folder of BED files. 
        Cannot be used together with -b flag.'''
    )

    # OPTIONAL: Lists all headers in a vcf then exits
    parser.add_argument(
        '-l', '--settings_list', action='store_true', 
        help='''Returns a list of all availabile headers from the VCF and the 
        source of the data, then exits.'''
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
    logger.info('running vcf_parse.py...\n{}'.format('---'*50))


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
        pt.load(args)
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
    logger.info('vcf_parse.py completed\n{}'.format('---'*50))
