#!/anaconda3/envs/python2/bin/python

"""
vcf_report.py

Object that deals with the loading and parsing of VCF files. 
Loaded as part of the vcf_parse.py program.

Author:     Erik Waskiewicz
Created:    31 Aug 2018
Edited:     Niamh Teague - 25 May 2023
Version:    0.1.1
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
            print(vcf_reader)
            vcf_records = []
            for var in vcf_reader:
                vcf_records.append(var)
            self.data = vcf_records
            self.logger.info('loading VCF completed')

        # load sample name from vcf
        self.sample = vcf_reader.samples[0]

        # load info, variant and format fields
        self.info_fields = vcf_reader.infos
        self.format_fields = vcf_reader.formats

        # load vep headers from vcf INFO field, split into list
        # added try except to fix issue #9
        try:
            self.vep_fields = self.info_fields['CSQ'][3].split(' ')[-1].split('|')
        except KeyError:
            self.vep_fields = [] # use empty list instead of None to avoid downstream errors

        # load output filepath
        if out is not None:
            self.output_dir = os.path.abspath(out)
        else:
            self.output_dir = os.path.abspath('.')
        self.report_path = os.path.join(
            self.output_dir, self.sample + '_VariantReport.txt')

        # make empty vep variable 
        self.config = None

    def load_config(self, config_file):
        """
        Load in config file that defines what annotations to include 
        within the variant report. Save config as a list.
        """
        self.logger.info('loading report config from {}'.format(
            os.path.abspath(config_file)))
        config = []
        with open(config_file, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            for line in reader:
                try:
                    config += [[ line[0], line[1], line[2] ]]
                except:
                    config += [[ line[0], line[1], '' ]]
            self.config = config
            self.logger.info('loading report config completed')


    def list_config(self):
        """
        Returns a list to screen containing all possible column headers
        and the source that the data comes from.
        """

        # Chrom, pos, alt, ref fields
        print('SYMBOL\tvep\tgene')
        print('CHROM\tvar\tchr')
        print('POS\tvar\tpos')
        print('REF\tvar\tref')
        print('ALT\tvar\talt')

        # preferred transcript and filter fields
        #print('Preferred\tpref')
        #print('Filter\tfilter')
        #print('Classification\tcustom')
        # format fields
        print('Frequency\tformat\tvaf')
        #for record in self.format_fields:
         #   print(record + '\tformat'

        # info fields - only want depth
        for record in self.info_fields:
            if record == 'DP':
                print('DP' + '\tinfo' + '\tdepth')

        # vep fields

        #print('dbSNP\tvep')
        #print('Cosmic\tvep')
        #print('HGMD\tvep')

        print('Feature\tvep\ttranscript')
        print('HGVSp\tvep\thgvs_p')
        print('HGVSc\tvep\thgvs_c')
        print('Consequence\tvep\tconsequence')
        print('EXON\tvep\texon')

 #       for record in self.vep_fields:
  #          if record == 'Feature'\
   #             or record == 'HGVSc'\
    #            or record == 'HGVSp'\
     #           or record == 'Consequence'\
      #          or record == 'EXON'\
       #         or record == 'SYMBOL':
        #        print(record + '\tvep' + '\tduh')


#    def make_variant_name(self, variant):
#        var_name = '{}{}{}>{}'.format(
 #           str(variant.CHROM),
  #          str(variant.POS),
   #         str(variant.REF),
    #        str(variant.ALT).strip('[]').replace(' ', '')
     #   )

#        return(var_name)

#    def parse_filter_field(self, variant):
 #       try:
  #          out = [str(variant.FILTER)]
   #     except:
    #        out = ['']
#
 #       return(out)


    def parse_chrom_field(self, variant):

        for variant in self.data:
            out = [variant.CHROM]
        else:
            out = ['']
            #print(out)
        return(out)


    def parse_pos_field(self, variant):
        for variant in self.data:
            out = [variant.POS]
        else:
            out = ['']

        return(out)


    def parse_ref_field(self, variant):
        for variant in self.data:
            out = [variant.REF]
        else:
            out = ['']

        return(out)


    def parse_alt_field(self, variant):
        for variant in self.data:
            out = variant.ALT
        else:
            out = ['']

        return(out)


    def parse_info_field(self, var, field):
        try:
            out = [str(var.INFO[field])]

        except:
            out = ['']

        return(out)


    def parse_format_field(self, variant, field):
        for sample in variant:
            if sample.sample == self.sample:
                try:
                    out = [ sample[field] ]
                except:
                    out = ['']
 
                # custom setting for allele freq
                if field == 'Frequency':
                    out = [ sample['AD'] ]
                    ref = float(out[0][0])
                    alt = float(out[0][1])
                    freq = float((alt / (ref + alt)) * 100)
                    out = ['{}%'.format(round(freq, 2))]

                # custom setting for genotype
                #if field == 'GT':
                 #   gt = out[0]
                  #  if gt == '0/1':
                   #     out = ['HET']
                    #if gt == '1/1':
                     #   out = ['HOM_VAR']
                    #if gt == '0/0':
                     #   out = ['HOM_REF']

        return(out)

    def parse_vep_field(self, field, vep):
        if vep:
            try:
                pos = self.vep_fields.index(field) #chrom?
                out = str(vep[pos]) #gene
            except:
                out = ''

            # custom exon/intron tweak
            if field in ('EXON', 'INTRON'):
                out = out.replace('/', '|')

            # custom HGVS coding/protein sequence tweak
            if field in ('HGVSc', 'HGVSp'):
                try:
                    split = out.split(':')
                    out = split[1]
                except:
                    pass

            # custom existing variantion field tweak
           # if field in ('dbSNP', 'Cosmic', 'HGMD'):
                # split existing variation field by & sign
            #    existing_variation_pos = self.vep_fields.index('Existing_variation')
             #   existing_variation = str(vep[existing_variation_pos]).split('&')

                # define identifier for each different annotation type
              #  if field == 'dbSNP':
               #     id = 'rs'
                #if field == 'Cosmic':
                 #   id = 'COSM'
                #if field == 'HGMD':
                 #   id = 'CM'

            # make output string containing only records of the desired type
#            out_list = ''
 #           for item in existing_variation:
#
 #               if item.startswith(id):
  #                  out_list += '{},'.format(str(item))
   #         out = out_list.rstrip(',')

            # custom exac/1kg tweak
            #if field in ('ExAC_AFR_MAF', 'ExAC_AMR_MAF', 'ExAC_EAS_MAF', 'ExAC_FIN_MAF',
             #   'ExAC_NFE_MAF', 'ExAC_SAS_MAF', 'ExAC_OTH_MAF', 'AFR_MAF', 'AMR_MAF',
              #  'EAS_MAF', 'EUR_MAF', 'SAS_MAF'):

               # out_string = ''

                #try:
                 #   for record in out.split('&'):
                  #      split = record.split(':')
                   #     percent = float(split[1]) * 100
                    #    out_string += '{}:{}%,'.format(split[0], str(percent))
                   # out = out_string.rstrip(',')
                #except:
                 #   pass

        else:
            out = 'No VEP output'

        return([out])


    def make_header(self):
        # gene to be first column
        # config file provided
        header = ''
        #if self.config:
        for annotation in self.config:
            if annotation[2] != '':
                header +=  '\t' + annotation[2]
            else:
                header +=  '\t' + annotation[0]

        # config file not provided - all headers
#        else:
 #           header += '\tPreferred\tClassification\tFilter'
  #          for annotation in self.info_fields:
   #             if annotation != 'CSQ':
    #                header += '\t' + annotation
     #       for annotation in self.format_fields:
      #          header += '\t' + annotation
       #     for annotation in self.vep_fields:
        #        header += '\t' + annotation

        # add newline and return
        header += '\n'
        return(header)


    def make_record_config(self, setting, variant, vep=None):
        """
        Makes a line of the variant report if config are present
        """
        out = ['']

        #  chrom
        if setting[2] == 'chr':
            out = self.parse_chrom_field(setting[0])

        #  pos
        if setting[2] == 'pos':
            out = self.parse_pos_field(setting[0])

        #  ref
        if setting[2] == 'ref':
            out = self.parse_ref_field(setting[0])

        #  alt
        if setting[2] == 'alt':
            out = self.parse_alt_field(setting[0])

        # preferred
        if setting[1] == 'pref':
            out = ['Unknown']

        # filter
        if setting[1] == 'filter':
            out = self.parse_filter_field(variant)

        # info
        if setting[1] == 'info':
            out = self.parse_info_field(variant, setting[0])

        # format
        if setting[1] == 'format':
            out = self.parse_format_field(variant, setting[0])

        # vep header
        if setting[1] == 'vep':
            out = self.parse_vep_field(setting[0], vep)

        return(out)


    def make_record_no_config(self, variant, vep=None):
        """
        Makes a line of the variant report if no config are present
        """
        out = []

        # preferred
        out += ['Unknown']
        out += ['']

        # filter
        out += self.parse_filter_field(variant)

        # info - don't include CSQ field, this is parsed as part of the vep parser
        for annotation in self.info_fields:
            if annotation != 'CSQ':
                out += self.parse_info_field(variant, annotation)

        # format
        for annotation in self.format_fields:
            out += self.parse_format_field(variant, annotation)

        # vep
        for annotation in self.vep_fields:
            out += self.parse_vep_field(annotation, vep)

        return(out)


    def make_report(self, filter_setting):
        """
        Contains a lot of nested loops, overview of loop structure:

        - open file to save output to
        - saves headers to output file
        - loops through each variant:
           - if variant has VEP annotation:
              - loop through each transcript:
                 - if config file provided: 
                    - loop through config and add to output list
                 - if no config: 
                    - loop through all annotations and add to output list
                 - save output list to output file
           - if no VEP annotations:
              - if config file provided: 
                 - loop through config and add to output list
              - if no config: 
                 - loop through all annotations and add to output list
              - save output list to output file
        - remove duplicate records  
        """
        self.logger.info('writing variant report')

        # open empty output file
        outfile = open(self.report_path, 'w')
        report_writer = csv.writer(outfile, delimiter='\t')

        # loop through variants
        for var in self.data:

            # PASS filter - pass will be empty - [], anything else will be filtered out
            if filter_setting and var.FILTER :
                pass

            else:
                # make variant name
                #variant = self.make_variant_name(var)

                # if VEP annotation exists, loop through each transcript
                try:
                    vep = var.INFO['CSQ']

                    for record in range(len(vep)):
                        out = []
                        vep_split = vep[record].split('|')

                        # filter out any transcripts that dont begin with NM
                        transcript_col = self.vep_fields.index('Feature')
                        if vep_split[transcript_col].startswith('NM'):

                            # if config file provided, parse annotations
                            if self.config:
                                for annotation in self.config:
                                    out += self.make_record_config(annotation, var, vep=vep_split)

                            # if no config file - include all annotations
                            # filter and preferred must be first, in that order
                            else:
                                out = self.make_record_no_config(var, vep=vep_split)

                            # save to file then repeat for all transcripts
                            report_writer.writerow(out)

                # if variant has no vep annotations
                except:
                    out = []

                    # if config file provided, parse annotations
                    if self.config:
                        for annotation in self.config:
                            out += self.make_record_config(annotation, var)

                    # if no config file - include all annotations
                    # filter and preferred must be first, in that order
                    else:
                        out = self.make_record_no_config(var)

                    # save to file then repeat for next variant
                    report_writer.writerow(out)

        # once loop has finished, close the output file
        outfile.close()

        # remove duplicates in output file
        uniq = os.popen('cat {} | uniq'.format(self.report_path)).read()

        # write final report
        out = open(self.report_path, 'w') 
        out.write(self.make_header())
        out.write(uniq)
        out.close()
        self.logger.info('variant report completed - {}'.format(self.report_path))
