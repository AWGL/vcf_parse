#!/anaconda3/envs/python2/bin/python

"""
vcf_report.py

Object that deals with the loading and intersecting of BED files. 
Loaded as part of the vcf_parse.py program.

Author:     Erik Waskiewicz
Created:    31 Aug 2018
Version:    1.0.0
Updated:    23 Oct 2018
"""


import os
import csv
import logging
import re


# set global variable to point to bedtools executables
BEDTOOLS_PATH = '/Users/erik/Applications/bedtools2/bin/'


# -- BED CLASS --------------------------------------------------------

class bed_object:
    def __init__(self):
        """
        Object properties that are loaded when the oject is created.
        The logger deals with all status messages from the object and
        is a child of the main logger vcf_parse.
        """
        self.logger = logging.getLogger('vcf_parse.bed')


    def make_intersect_bed(self, bedfile, in_vcf):
        """
        - Takes a variant report and turns into a BED file by parsing 
          the variant description, and saves file as <report>.temp
        - Intersects the temp and input BED files
        - Outputs intersect BED containing variants from the report 
          that fall within the BED file. This file also contains the 
          original variant description so that the variant report can 
          easily be filtered in the next function.
        """
        # open variant report to loop through
        with open(in_vcf.report_path, 'r') as report:
            reader = csv.reader(report, delimiter='\t')
            # open empty temp file to write output to
            with open(in_vcf.report_path + '.temp', 'w') as out: 
                # for each line, split variant into BED format and save
                for line in reader:
                    if line[0][0] != '#':
                        variant = line[1].split(':')
                        pos = re.sub('[^0-9]', '', variant[1])
                        out.write('{}\t{}\t{}\t{}\n'.format(
                            variant[0], str(int(pos) - 1), pos, line[1]
                        ))
            out.close()
        report.close()

        # make output folder if it doesnt exist, based on input folder name
        in_folder = os.path.split(os.path.dirname(bedfile))[-1]
        out_folder = os.path.join(in_vcf.output_dir, in_folder)
        if not os.path.exists(out_folder):
            os.mkdir(out_folder)

        # intersect report bed and input bed, save as variable
        report_bed = in_vcf.report_path + '.temp'
        intersect_command = (
            '{}intersectBed -a {} -b {}'.format(BEDTOOLS_PATH, report_bed, bedfile)
        )
        try:
            results_intersect = os.popen(intersect_command).read()
        except IOError:
            self.logger.warn('BEDTools error')

        # write bed variable to file, remove report bed
        self.bed_name = os.path.basename(bedfile).rstrip('.bed')
        self.intersect_bed = os.path.join(
            out_folder, '{}_{}_intersect.bed'.format(
                in_vcf.sample, self.bed_name))
        self.out_folder = out_folder

        out = open(self.intersect_bed, 'w') 
        out.write(results_intersect)
        out.close()
        os.remove(report_bed)


    def apply_bed(self, bedfile, in_vcf):
        """
        Takes an intersect BED from the make_intersect_bed function,
        compares the variant ID within this to the variant ID within 
        the variant report, if the two match, keeps the line of the 
        report, otherwise discards it
        """
        # load temp intersect bed, read variants into list
        with open(self.intersect_bed, 'r') as bed:
            results = csv.reader(bed, delimiter='\t')
            keep = []
            for line in results:
                keep.append(line[3])

        # open empty file
        outfile = os.path.join(
            self.out_folder, '{}_{}_VariantReport.txt'.format(
                in_vcf.sample, self.bed_name))
        bed = open(outfile, 'w')
        bed_report = csv.writer(bed, delimiter='\t')

        # loops through original report, keeps if there is a match with
        # the bed list. Variant description from the report was kept in
        # the intersect BED, so they can be directly compared
        with open(in_vcf.report_path) as report:
            results = csv.reader(report, delimiter='\t')
            for line in results:
                if line[0][0] == '#':
                    bed_report.writerow(line)
                if line[1] in keep:
                    bed_report.writerow(line)
        
        # log and remove intersect BED
        self.logger.info('applied BED file - {}'.format(outfile))
        os.remove(os.path.abspath(self.intersect_bed))


    def apply_single(self, bedfile, in_vcf):
        """
        Takes a single BED file and calls the functions required to 
        apply it
        """
        self.logger.info('loading BED file 1 of 1: {}'.format(
            os.path.abspath(bedfile)))
        self.make_intersect_bed(bedfile, in_vcf)
        self.apply_bed(bedfile, in_vcf)

    
    def apply_multiple(self, bed_folder, in_vcf):
        """
        Takes a folder of BED files, loops through each one in turn and
        calls the functions required to apply it
        """
        n = len(os.listdir(bed_folder))
        for i, file in enumerate(os.listdir(bed_folder)):
            in_bed = os.path.join(bed_folder, file)
            self.logger.info('loading BED file {} of {}: {}'.format(
                i+1, n, os.path.abspath(in_bed)))
            self.make_intersect_bed(in_bed, in_vcf)
            self.apply_bed(in_bed, in_vcf)
