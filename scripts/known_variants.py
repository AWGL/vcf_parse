#!/anaconda3/envs/python2/bin/python

"""
known_variants.py

Object that deals with loading and annotating known variants. 
Loaded as part of the vcf_parse.py program.

Author:     Erik Waskiewicz
Created:    30 Oct 2018
Version:    0.1.0
Updated:    30 Oct 2018
"""


import os
import csv
import vcf
import logging


# -- KNOWN VARIANTS CLASS ---------------------------------------------

class known_variants:
    def __init__(self):
        """
        Object properties that are loaded when the oject is created.
        The logger deals with all status messages from the object and
        is a child of the main logger vcf_parse.
        """
        self.logger = logging.getLogger('vcf_parse.known')


    def load_known_variants(self, inp):
        """
        Load in vcf and save as list
        """
        # read input vcf with pyvcf package, save as list
        self.logger.info(
            'loading known variants from {}'.format(os.path.abspath(inp)))

        with open(inp, 'r') as vcf_input:
            vcf_reader = vcf.Reader(vcf_input)
            vcf_list = []
            vcf_records = []
            for var in vcf_reader:
                var_name = '{}:{}{}>{}'.format(
                    str(var.CHROM), 
                    str(var.POS), 
                    str(var.REF), 
                    str(var.ALT).strip('[]').replace(' ', '')
                )
                classification = var.INFO['Classification']
                vcf_list.append(var_name)
                vcf_records.append([var_name, classification])

            self.list = vcf_list
            self.records = vcf_records
            self.logger.info('loading known variants completed')


    def apply_known_variants(self, report):
        """
        check that there is a Classification column
        Compares variant id with list of known variants, annotates if there is a match
        """
        # set report path
        report_path = report.report_path

        # find classification column in config file
        if report.annotations:
            for record in report.annotations:
                if record[0] == 'Classification':
                    if record[2] != '':
                        classification_id = record[2]
                    else:
                        classification_id = record[0]
        else:
            classification_id = 'Classification'

        if self.list:
            # open report file and new temp file to save output
            report_temp = os.path.join(report_path + '.temp')
            f1 = open(report_path, 'rb')
            reader = csv.reader(f1, delimiter='\t')
            f2 = open(report_temp, 'wb')
            writer = csv.writer(f2, delimiter='\t')

            # load header, 
            header = next(reader)
            variant_column = 1

            # find classification column, make one if not present
            try:
                classification_column = header.index(classification_id)
            except:
                classification_column = len(header)
                header += ['Classification']

            writer.writerow(header)

            # loop though the report, annotate with classification if there is a match
            for row in reader:
                if row[variant_column] in self.list:
                    i = ''
                    for r in self.records:
                        if r[0] == row[variant_column]:
                            i += '{},'.format(r[1])
                    writer.writerow(row[0:classification_column] + [i.rstrip(',')] + row[classification_column+1:])
                else:
                    writer.writerow(row)

            # tidy up
            self.logger.info('known variants applied')
            f1.close()
            f2.close()
            os.remove(report_path)
            os.rename(report_temp, report_path)
        
        # if error, skip adding known variants
        else:
            self.logger.warn(
                'could not load known variants file provided, skipping step.'
            )
            return
