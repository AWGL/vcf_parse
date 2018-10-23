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
        '''
        Object properties that are loaded when the oject is created.
        The logger deals with all status messages from the object.
        '''
        self.logger = logging.getLogger('vcf_parse.vcf')


    def load_data(self, inp, out):
        # read input vcf with pyvcf package, save as list
        self.logger.info('loading VCF file from {}'.format(os.path.abspath(inp)))
        with open(inp, 'r') as vcf_input:
            vcf_reader = vcf.Reader(vcf_input)
            vcf_records = []
            for var in vcf_reader:
                vcf_records.append(var)
            self.data = vcf_records
            self.logger.info('loading VCF completed')

        # load sample name from vcf
        self.sample = vcf_reader.samples[0]

        # load vep headers from vcf INFO field, split into list
        self.info_fields = vcf_reader.infos
        self.vep_fields = self.info_fields['CSQ'][3].split(' ')[-1].split('|')
        self.format_fields = vcf_reader.formats

        # load output filepath
        if out is not None:
            self.output_dir = os.path.abspath(out)
        else:
            self.output_dir = os.path.abspath('.')
        self.report_path = os.path.join(self.output_dir, self.sample + '_VariantReport.txt')

        # make empty vep variable 
        self.annotations = None


    def settings(self, settings_file):
        # Load in external file that defines what vep annotations to use
        self.logger.info('loading report settings from {}'.format(os.path.abspath(settings_file)))
        settings = []
        with open(settings_file, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            for line in reader:
                settings += [[ line[0], line[1] ]]
            self.annotations = settings
            self.logger.info('loading report settings completed')


    def list_settings(self):
        '''
        Returns a list to screen containing all possible column headers
        and the source that the data comes from.
        '''
        # info fields
        for record in self.info_fields:
            print(record + '\tinfo')
        # format fields
        for record in self.format_fields:
            print(record + '\tformat')
        # preferred transcript and filter fields
        print('Preferred\tpref')
        print('Filter\tfilter')
        # vep fields
        for record in self.vep_fields:
            print(record + '\tvep')


    def make_report(self):
        self.logger.info('writing variant report')

        # open empty output file
        outfile1 = open(self.report_path, 'w')
        outfile = csv.writer(outfile1, delimiter='\t')

        # ----- report body --------
        for var in self.data:
            # make variant name
            variant = str(var.CHROM) + ':' + str(var.POS) + str(var.REF) + '>' + str(var.ALT).strip('[]').replace(' ', '')
            
            
            # -----------------------------------------
            # test if vep annotations present
            try:
                # If VEP annotation exists, loop through each transcript
                vep = var.INFO['CSQ']
                # -----------------------------------------
                # if annotations list included, loops through the list and include those annotations
                for record in range(len(vep)):
                    out = []
                    vep_split = vep[record].split('|')

                    # -----------------------------------------
                    # own annotation list
                    if self.annotations:
                        for annotation in self.annotations:
                            if annotation[1] == 'pref':
                                out += ['Unknown']
                            if annotation[1] == 'filter':
                                try:
                                    out += [str(var.FILTER)]
                                except:
                                    out += ['']
                            if annotation[1] == 'info':
                                try:
                                    out += [str(var.INFO[annotation[0]])]
                                except:
                                    out += ['']
                            if annotation[1] == 'format':
                                for sample in var:
                                    if sample.sample == self.sample:
                                        try:
                                            out += [ sample[ annotation[0] ] ]
                                        except:
                                            out += ['']
                            if annotation[1] == 'vep':
                                try:
                                    pos = self.vep_fields.index(annotation[0])
                                    out += [str(vep_split[pos])]
                                except:
                                    out += ['']

                    # -----------------------------------------
                    # no annotation list
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
                    
                    # -----------------------------------------
                    # write to file
                    outfile.writerow([self.sample] + [variant] + out)

            # -----------------------------------------
            # no vep_output 
            except:
                out = []

                # -----------------------------------------
                # own annotation list
                if self.annotations:
                    for annotation in self.annotations:
                        if annotation[1] == 'pref':
                            out += ['Unknown']
                        if annotation[1] == 'filter':
                            try:
                                out += [str(var.FILTER)]
                            except:
                                out += ['']
                        if annotation[1] == 'info':
                            try:
                                out += [str(var.INFO[annotation[0]])]
                            except:
                                out += ['']
                        if annotation[1] == 'format':
                            for sample in var:
                                if sample.sample == self.sample:
                                    try:
                                        out += [ sample[ annotation[0] ] ]
                                    except:
                                        out += ['']
                        if annotation[1] == 'vep':
                            out += ['No VEP output']
                
                # -----------------------------------------
                # all annotations
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
                        out += ['No VEP output']

                # -----------------------------------------
                # write to file
                outfile.writerow([self.sample] + [variant] + out)

        outfile1.close()

        # ----- report header --------
        # Sample and variant are always the first two columns
        header = '#Sample\tVariant'

        # Customisible headers
        if self.annotations:
            for annotation in self.annotations:
                header +=  '\t' + annotation[0]

        # all headers
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
        
        # remove duplicates
        uniq = os.popen('cat ' + self.report_path + ' | sort | uniq').read()

        out = open(self.report_path, 'w') 

        out.write(header)
        out.write(uniq)
        out.close()
        self.logger.info('variant report completed - {}'.format(self.report_path))

