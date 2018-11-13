#!/anaconda3/envs/python2/bin/python

"""
bed_object.py

Object that deals with the loading and intersecting of BED files. 
Loaded as part of the vcf_parse.py program.

Author:     Erik Waskiewicz
Created:    31 Aug 2018
Version:    0.1.0
Updated:    31 Oct 2018
"""


import os
import csv
import logging
import re


# set global variable to point to bedtools executables - if using conda leave as empty string
BEDTOOLS_PATH = ''


# -- BED CLASS --------------------------------------------------------

class bed_object:
    def __init__(self):
        """
        Object properties that are loaded when the oject is created.
        The logger deals with all status messages from the object and
        is a child of the main logger vcf_parse.
        """
        self.logger = logging.getLogger('vcf_parse.bed')


    def make_report_bed(self, in_vcf):
        """
        Takes a variant report and turns into a BED file by parsing 
        the variant description, and saves file as <report>.temp
        """
        # open variant report to loop through
        with open(in_vcf.report_path, 'r') as report:
            reader = csv.reader(report, delimiter='\t')

            # open empty temp file to write output to
            with open(in_vcf.report_path + '.temp', 'w') as out: 

                # for each line, split variant into BED format and save
                for line in reader:

                    if line[0] != 'SampleID':


                        variant = line[1].split(':')
                        ref = variant[1].split('>')[0].strip('0123456789')
                        alt = variant[1].split('>')[1]

                        #Account for indels overlapping gene bed
                        if len(ref) > 1:

                            start_pos = int(variant[1].strip('AGTC>')) - 1
                            end_pos = start_pos + len(ref) + 1

                        else:
                
                            start_pos = int(variant[1].strip('AGTC>')) - 1
                            end_pos =  variant[1].strip('AGTC>')


                        out.write('{}\t{}\t{}\t{}\n'.format(
                            variant[0],
                            start_pos,
                            end_pos,
                            line[1]
                        ))
            out.close()
        report.close()


    def make_intersect_bed(self, bedfile, in_vcf):
        """
        - Intersects the report BED and input BED files
        - Outputs intersect BED containing variants from the report 
          that fall within the BED file. This file also contains the 
          original variant description so that the variant report can 
          easily be filtered in the next function.
        """
        # intersect report bed and input bed, save as variable
        report_bed = in_vcf.report_path + '.temp'
        intersect_command = (
            '{}intersectBed -a {} -b {}'.format(BEDTOOLS_PATH, report_bed, bedfile)
        )
        try:
            results_intersect = os.popen(intersect_command).read()
        except IOError:
            self.logger.warn('BEDTools error')

        # write bed variable to file
        self.bed_name = os.path.basename(bedfile).split('.')[0]
        self.intersect_bed = os.path.join(
            in_vcf.output_dir, '{}_{}_intersect.bed'.format(
                in_vcf.sample, self.bed_name))
        #self.out_folder = out_folder

        out = open(self.intersect_bed, 'w') 
        out.write(results_intersect)
        out.close()


    def apply_bed(self, bedfile, in_vcf, out_folder):
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
            out_folder, '{}_{}_VariantReport.txt'.format(
                in_vcf.sample, self.bed_name))
        bed = open(outfile, 'w')
        bed_report = csv.writer(bed, delimiter='\t')

        # loops through original report, keeps if there is a match with
        # the bed list. Variant description from the report was kept in
        # the intersect BED, so they can be directly compared
        with open(in_vcf.report_path) as report:
            results = csv.reader(report, delimiter='\t')
            for line in results:
                if line[0] == 'SampleID':
                    bed_report.writerow(line)
                if line[1] in keep:
                    bed_report.writerow(line)
        
        # log and remove intersect BED
        self.logger.info('applied BED file - {}'.format(outfile))
        os.remove(os.path.abspath(self.intersect_bed))


    def apply_single(self, bedfile, in_vcf):
        """
        Takes a single BED file and calls the functions required to 
        apply it. Saves to the same directory as the variant report.
        """
        self.logger.info('loading BED file 1 of 1: {}'.format(
            os.path.abspath(bedfile)))
        
        # make temporary report BED file
        self.make_report_bed(in_vcf)

        # make intersect BED and apply to variant report
        self.make_intersect_bed(bedfile, in_vcf)
        self.apply_bed(bedfile, in_vcf, in_vcf.output_dir)

        # remove temp report BED
        os.remove(in_vcf.report_path + '.temp')

    
    def apply_multiple(self, bed_folder, in_vcf):
        """
        Takes a folder of BED files, loops through each one in turn and
        calls the functions required to apply it.
        Makes a new directory within the output directory, with the 
        same name as the input BED folder, to save the output
        """
        # make output folder if it doesnt exist, based on input folder name
        in_folder = os.path.basename(bed_folder)
        out_folder = os.path.join(in_vcf.output_dir, in_folder)
        if not os.path.exists(out_folder):
            os.mkdir(out_folder)

        # make temporary report BED file once
        self.make_report_bed(in_vcf)

        # total number of files for logger
        n = len([name for name in os.listdir(bed_folder) 
            if os.path.isfile(os.path.join(bed_folder, name))])
        
        # loop through BED files within folder
        i = 1
        for file in os.listdir(bed_folder):
            in_bed = os.path.join(bed_folder, file)
            if os.path.isfile(in_bed):
                self.logger.info('loading BED file {} of {}: {}'.format(
                    i, n, os.path.abspath(in_bed)))

                # make intersect BED and apply to variant report
                self.make_intersect_bed(in_bed, in_vcf)
                self.apply_bed(in_bed, in_vcf, out_folder)

                i+=1

        # remove temp report BED
        os.remove(in_vcf.report_path + '.temp')
