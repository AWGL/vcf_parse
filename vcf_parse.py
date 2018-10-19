import argparse

from scripts.vcf_report import vcf_report
from scripts.preferred_transcripts import preferred_transcripts
from scripts.bed_object import bed_object


## -- PARSE INPUT ARGUMENTS -------------------------------------------

def get_args():
    '''
    Use argparse package to take arguments from the command line and 
    store them as an argparse object.
    '''

    parser = argparse.ArgumentParser(
        description='''Takes a VCF file and converts it into a tsv format 
        variant report. Can also take a BED file or a folder of BED 
        files and output a seperate variant report containing only the 
        calls within each of the BED files.'''
    )

    # required - vcf file
    parser.add_argument(
        'input', action='store', 
        help='Filepath to input vcf file. Required.'
    )

    # optional - output folder, defaults to current directory if empty
    parser.add_argument(
        '-O', '--output', action='store', 
        help='''Filepath to folder where output reports will be saved. If 
        empty, defaults to current directory.'''
    )

    # optional - file containing the headers for the report
    parser.add_argument(
        '-s', '--settings', action='store', 
        help='''Filepath to settings file. This is a tab seperated file which 
        specifies the order of the report headers in the first column and the 
        location to find the data in the second column. Column headers must 
        match up with how they appear in the VCF, for a full list of the 
        available settings, run vcf_parse with the -l flag.'''
    )

    # optional - list of preferred transcripts
    parser.add_argument(
        '-t', '--transcripts', action='store', 
        help='''Filepath to preferred transcripts file. must be a tab 
        seperated file with preferred transcripts in the seond column. If 
        empty, all entries in the preferred transcript column (if present) 
        will be false.'''
    )
    
    # optional -  either a single bed file or a folder containing bed 
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

    # optional - lists the VEP headers in a vcf, if this flag is used 
    # then it will return the list and exit
    parser.add_argument(
        '-l', '--settings_list', action='store_true', 
        help='Returns a list of VEP headers from the VCF then exits.'
    )

    return parser.parse_args()


# -- CALL FUNCTIONS ---------------------------------------------------

if __name__ == '__main__':
    args = get_args()
    report = vcf_report()

    # Load data using load_data function within vcf_report.py
    report.load_data(args.input, args.output)

    # print vep headers and exit if -l flag called
    if args.settings_list:
        for record in report.vep_fields:
            print(record + '\tvep')
        exit()
    
    # get vep settings
    if args.settings:
        report.settings(args.settings)
    else:
        print('No settings file found -- outputting all data from VCF.')

    # make total report
    report.make_report()

    # apply preferred transcripts
    if args.transcripts:
        pt = preferred_transcripts()
        pt.load(args)
        pt.apply(report.report_path)
    else:
        print('''No preferred transcripts found -- preferred transcripts 
            column will all be labelled as FALSE.''')

    # make single bed report
    if args.bed:
        bed = bed_object()
        bed.apply_single(args.bed, report)

    # make multiple bed report
    if args.bed_folder:
        bed = bed_object()
        bed.apply_multiple(args.bed_folder, report)

