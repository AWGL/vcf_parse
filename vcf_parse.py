import argparse

from scripts.vcf_report import vcf_report
from scripts.preferred_transcripts import preferred_transcripts
from scripts.bed_object import bed_object


# ----------------- PARSE ARGUMENTS -----------------------------------

# use argparse package to take arguments from the command line
def get_args():
    parser = argparse.ArgumentParser()

    # standard inputs
    parser.add_argument('-I', '--input', action='store', required=True, help='Filepath to input vcf file. Required.')
    parser.add_argument('-O', '--output', action='store', help='Filepath to output location. Defaults to current directory.')
    parser.add_argument('-v', '--vep_settings', action='store', help='Filepath to vep_settings file.')
    parser.add_argument('-l', '--vep_list', action='store_true', help='Returns a list of VEP headers from the VCF then exits.')
    parser.add_argument('-t', '--transcripts', action='store', help='Filepath to preferred transcripts file.')
    
    # only one of these can be used
    bed_files = parser.add_mutually_exclusive_group()
    bed_files.add_argument('-b', '--bed', action='store', help='Filepath to bed file.')
    bed_files.add_argument('-B', '--bed_folder', action='store', help='Filepath to folder of bed files.')

    # not used yet
    parser.add_argument('-m', '--multisample', action='store_true', help='VCF file contains multiple samples.')

    return parser.parse_args()


# ----------------- CALL FUNCTIONS ------------------------------------

if __name__ == '__main__':
    args = get_args()
    report = vcf_report()

    # get data
    report.load_data(args.input, args.output)

    # print vep headers and exit if -l flag called
    if args.vep_list:
        for record in report.vep_fields:
            print(record)
        exit()
    
    # get vep settings
    if args.vep_settings:
        report.vep_settings(args.vep_settings)

    # make total report
    report.make_report()

    # apply preferred transcripts
    if args.transcripts:
        pt = preferred_transcripts()
        pt.load(args)
        pt.apply(report.report_path)
    else:
        print('No preferred transcripts found')

    # make single bed report
    if args.bed:
        bed = bed_object()
        bed.apply_single(args.bed, report)

    # make multiple bed report
    if args.bed_folder:
        bed = bed_object()
        bed.apply_multiple(args.bed_folder, report)

