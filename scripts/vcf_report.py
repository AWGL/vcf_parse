#!/anaconda3/envs/python2/bin/python

"""
vcf_report.py

Object that deals with the loading and parsing of VCF files. 
Loaded as part of the vcf_parse.py program.

Author:     Erik Waskiewicz
Created:    31 Aug 2018
Version:    1.0.0
Updated:    23 Oct 2018
"""


import os
import vcf
import csv
import logging


# ----------------- REPORT CLASS --------------------------------------
class vcf_report:
    def __init__(self):
        """
        Object properties that are loaded when the oject is created.
        The logger deals with all status messages from the object and
        is a child of the main logger vcf_parse.
        """
        self.logger = logging.getLogger('vcf_parse.vcf')


    def load_data(self, inp, out):
        """Load in data from a VCF"""
        # read input vcf with pyvcf package, save as list
        self.logger.info(
            'loading VCF file from {}'.format(os.path.abspath(inp)))
        with open(inp, 'r') as vcf_input:
            vcf_reader = vcf.Reader(vcf_input)
            vcf_records = []
            for var in vcf_reader:
                vcf_records.append(var)
            self.data = vcf_records
            self.logger.info('loading VCF completed')

        # load sample name from vcf
        self.sample = vcf_reader.samples[0]

        # load info and format fields
        self.info_fields = vcf_reader.infos
        self.format_fields = vcf_reader.formats

        # load vep headers from vcf INFO field, split into list
        self.vep_fields = self.info_fields['CSQ'][3].split(' ')[-1].split('|')

        # load output filepath
        if out is not None:
            self.output_dir = os.path.abspath(out)
        else:
            self.output_dir = os.path.abspath('.')
        self.report_path = os.path.join(
            self.output_dir, self.sample + '_VariantReport.txt')

        # make empty vep variable 
        self.annotations = None


    def settings(self, settings_file):
        """
        Load in settings file that defines what annotations to include 
        within the variant report. Save settings as a list.
        """
        self.logger.info('loading report settings from {}'.format(
            os.path.abspath(settings_file)))
        settings = []
        with open(settings_file, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            for line in reader:
                settings += [[ line[0], line[1] ]]
            self.annotations = settings
            self.logger.info('loading report settings completed')


    def list_settings(self):
        """
        Returns a list to screen containing all possible column headers
        and the source that the data comes from.
        """
        # preferred transcript and filter fields
        print('Preferred\tpref')
        print('Filter\tfilter')

        # info fields
        for record in self.info_fields:
            print(record + '\tinfo')

        # format fields
        for record in self.format_fields:
            print(record + '\tformat')

        # vep fields
        for record in self.vep_fields:
            print(record + '\tvep')


    def make_report(self):
        """
        Contains a lot of nested loops, overview of loop structure:

        - open file to save output to
        - saves headers to output file
        - loops through each variant:
           - if variant has VEP annotation:
              - loop through each transcript:
                 - if settings file provided: 
                    - loop through settings and add to output list
                 - if no settings: 
                    - loop through all annotations and add to output list
                 - save output list to output file
           - if no VEP annotations:
              - if settings file provided: 
                 - loop through settings and add to output list
              - if no settings: 
                 - loop through all annotations and add to output list
              - save output list to output file
        - remove duplicate records  
        """
        self.logger.info('writing variant report')

        # -- report header --------------------------------------------
        # Sample and variant are always the first two columns
        header = '#Sample\tVariant'

        # settings file provided
        if self.annotations:
            for annotation in self.annotations:
                header +=  '\t' + annotation[0]

        # settings file not provided - all headers
        else:
            header += '\tFilter\tPreferred'
            for annotation in self.info_fields:
                if annotation != 'CSQ':
                    header += '\t' + annotation
            for annotation in self.format_fields:
                header += '\t' + annotation
            for annotation in self.vep_fields:
                header += '\t' + annotation
        header += '\n'

        # -- report body ----------------------------------------------
        # open empty output file
        outfile = open(self.report_path, 'w')
        report_writer = csv.writer(outfile, delimiter='\t')

        # loop through variants
        for var in self.data:

            # make variant name
            variant = '{}:{}{}>{}'.format(
                str(var.CHROM), 
                str(var.POS), 
                str(var.REF), 
                str(var.ALT).strip('[]').replace(' ', '')
            )
            
            # if VEP annotation exists, loop through each transcript
            try:
                vep = var.INFO['CSQ']
                for record in range(len(vep)):
                    out = []
                    vep_split = vep[record].split('|')

                    # if settings file provided, parse annotations
                    if self.annotations:
                        for annotation in self.annotations:

                            # preferred
                            if annotation[1] == 'pref':
                                out += ['Unknown']

                            # filter
                            if annotation[1] == 'filter':
                                try:
                                    out += [str(var.FILTER)]
                                except:
                                    out += ['']

                            # info
                            if annotation[1] == 'info':
                                try:
                                    out += [str(var.INFO[annotation[0]])]
                                except:
                                    out += ['']

                            # format
                            if annotation[1] == 'format':
                                for sample in var:
                                    if sample.sample == self.sample:
                                        try:
                                            out += [ sample[annotation[0]] ]
                                        except:
                                            out += ['']

                            # vep header
                            if annotation[1] == 'vep':
                                try:
                                    pos = self.vep_fields.index(annotation[0])
                                    out += [str(vep_split[pos])]
                                except:
                                    out += ['']

                    # if no settings file - include all annotations
                    # filter and preferred must be first, in that order
                    else:
                        # filter
                        try:
                            out += [var.FILTER]
                        except:
                            out += ['']

                        # preferred
                        out += ['Unknown']

                        # info
                        for annotation in self.info_fields:
                            if annotation != 'CSQ':
                                try:
                                    out += [var.INFO[annotation]]
                                except:
                                    out += ['']

                        # format
                        for annotation in self.format_fields:
                            for sample in var:
                                if sample.sample == self.sample:
                                    try:
                                        out += [ sample[ annotation ] ]
                                    except:
                                        out += ['']

                        # vep
                        for annotation in self.vep_fields:
                            try:
                                pos = self.vep_fields.index(annotation)
                                out += [str(vep_split[pos])]
                            except:
                                out += ['']
                    
                    # save to file then repeat for all transcripts
                    report_writer.writerow([self.sample] + [variant] + out)

            # if variant has no vep annotations
            except:
                out = []

                # if settings file provided, parse annotations
                if self.annotations:
                    for annotation in self.annotations:
                        # preferred
                        if annotation[1] == 'pref':
                            out += ['Unknown']

                        # filter
                        if annotation[1] == 'filter':
                            try:
                                out += [str(var.FILTER)]
                            except:
                                out += ['']

                        # info
                        if annotation[1] == 'info':
                            try:
                                out += [str(var.INFO[annotation[0]])]
                            except:
                                out += ['']

                        # format
                        if annotation[1] == 'format':
                            for sample in var:
                                if sample.sample == self.sample:
                                    try:
                                        out += [ sample[annotation[0]] ]
                                    except:
                                        out += ['']

                        # vep header
                        if annotation[1] == 'vep':
                            out += ['No VEP output']
                
                # if no settings file - include all annotations
                # filter and preferred must be first, in that order
                else:
                    # filter
                    try:
                        out += [var.FILTER]
                    except:
                        out += ['']

                    # preferred
                    out += ['Unknown']

                    # info
                    for annotation in self.info_fields:
                        if annotation != 'CSQ':
                            try:
                                out += [var.INFO[annotation]]
                            except:
                                out += ['']

                    # format
                    for annotation in self.format_fields:
                        for sample in var:
                            if sample.sample == self.sample:
                                try:
                                    out += [sample[annotation]]
                                except:
                                    out += ['']

                    # vep
                    for annotation in self.vep_fields:
                        out += ['No VEP output']

                # save to file then repeat for next variant
                report_writer.writerow([self.sample] + [variant] + out)

        # once loop has finished, close the output file
        outfile.close()
        
        # remove duplicates in output file
        uniq = os.popen('cat ' + self.report_path + ' | sort | uniq').read()

        # write final report
        out = open(self.report_path, 'w') 
        out.write(header)
        out.write(uniq)
        out.close()
        self.logger.info('variant report completed - {}'.format(self.report_path))
