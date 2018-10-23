#!/anaconda3/envs/python2/bin/python

"""
test.py

Unit tests for the vcf_parse.py program.

Author:     Erik Waskiewicz
Created:    22 Oct 2018
Version:    1.0.0
Updated:    23 Oct 2018
"""


import unittest
import os
import csv

from scripts.vcf_report import vcf_report
from scripts.preferred_transcripts import preferred_transcripts
from scripts.bed_object import bed_object


class TestVCF(unittest.TestCase):
    def setUp(self):
        """load in common files"""
        self.report = vcf_report()
        self.report.load_data(
            os.path.abspath('test/test.vcf'), os.path.abspath('test/')
            )
        self.report.make_report()


    def tearDown(self):
        """remove output files after test has run"""
        os.remove(self.report.report_path)
        self.report = None
        self.pt = None
        if os.path.isfile('test/test_bed_files/SAMPLE1_bed1_VariantReport.txt'):
            os.remove(os.path.abspath(
                'test/test_bed_files/SAMPLE1_bed1_VariantReport.txt'
                ))
        if os.path.isfile('test/test_bed_files/SAMPLE1_bed2_VariantReport.txt'):
            os.remove(os.path.abspath(
                'test/test_bed_files/SAMPLE1_bed2_VariantReport.txt'
                ))


    def test_vcf_parser_number_variants(self):
        """
        Check that number of variants loaded from VCF is correct, 
        should be 96
        """
        self.assertEqual(len(self.report.data), 96, 
            'Number of variants incorrect'
            )


    def test_vcf_parser_sample_name(self):
        """
        Check that sample name loaded from VCF is correct, 
        should be SAMPLE1
        """
        self.assertEqual(self.report.sample, 'SAMPLE1', 
            'Sample name incorrect'
            )


    def test_vcf_parser_info_fields(self):
        """
        Check that info fields loaded from VCF are correct, 
        should be the same as the expected list below
        """
        expected_info_fields = ['AC', 'AF', 'AN', 'BaseQRankSum',
        'DP', 'EXON', 'FC', 'GI', 'LCRLen', 'MQ', 'MQ0',
        'MQRankSum', 'TI', 'CSQ']

        report_info_fields = []
        for field in self.report.info_fields:
            report_info_fields.append(field)
        self.assertEqual(report_info_fields, expected_info_fields,
            'info fields incorrect'
            )


    def test_vcf_parser_format_fields(self):
        """
        Check that format fields loaded from VCF are correct, 
        should be the same as the expected list below
        """
        expected_format_fields = ['AD', 'GQ', 'GT', 'NL', 'SB', 'VF']
        report_format_fields = []
        for field in self.report.format_fields:
            report_format_fields.append(field)
        self.assertEqual(report_format_fields, expected_format_fields, 
            'format fields incorrect'
            )


    def test_vcf_parser_vep_fields(self):
        """
        Check that vep fields loaded from VCF are correct, 
        should be the same as the expected list below
        """
        expected_vep_fields = ['Allele', 'Consequence', 'IMPACT', 'SYMBOL', 
        'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 
        'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position',
        'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM',
        'DISTANCE', 'STRAND', 'FLAGS', 'VARIANT_CLASS', 'SYMBOL_SOURCE',
        'HGNC_ID', 'CANONICAL', 'TSL', 'APPRIS', 'CCDS', 'ENSP', 'SWISSPROT',
        'TREMBL', 'UNIPARC', 'REFSEQ_MATCH', 'GENE_PHENO', 'SIFT', 'PolyPhen',
        'DOMAINS', 'HGVS_OFFSET', 'GMAF', 'AFR_MAF', 'AMR_MAF', 'EAS_MAF', 
        'EUR_MAF', 'SAS_MAF', 'AA_MAF', 'EA_MAF', 'ExAC_MAF', 'ExAC_Adj_MAF',
        'ExAC_AFR_MAF', 'ExAC_AMR_MAF', 'ExAC_EAS_MAF', 'ExAC_FIN_MAF', 
        'ExAC_NFE_MAF', 'ExAC_OTH_MAF', 'ExAC_SAS_MAF', 'CLIN_SIG', 'SOMATIC',
        'PHENO', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 
        'MOTIF_SCORE_CHANGE']

        self.assertEqual(self.report.vep_fields, expected_vep_fields)


    def test_variant_report_number_variants(self):
        """
        Check that number of rows in the variant report is correct,
        should be 477
        """
        report_sum = sum(1 for line in open(os.path.abspath(
            'test/SAMPLE1_VariantReport.txt'
            )))
        self.assertEqual(report_sum, 477, 
            'Number of variants incorrect'
            )


    def test_load_settings(self):
        """
        Check that settings file is loaded correctly, the headers loaded in 
        should be the same as the expected list below
        """
        expected_settings = ['#Sample', 'Variant', 'Filter', 'GT', 'AF',
        'DP', 'Feature', 'Preferred', 'CANONICAL', 'Allele', 
        'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'DP', 'Feature_type']

        self.report.settings(os.path.abspath('test/config.txt'))
        self.report.make_report()
        with open(self.report.report_path) as f:
            reader = csv.reader(f, delimiter='\t')
            self.assertEqual(reader.next(), expected_settings)


    def test_preferred_transcripts_high_strictness_true(self):
        """
        Check that preferred transcripts are labelled correctly. 
        NM_001007553 should be true with high and low strictness
        """
        self.pt = preferred_transcripts()
        self.pt.load(os.path.abspath('test/PreferredTranscripts.txt'))
        self.pt.apply(self.report.report_path, 'high')
        with open(self.report.report_path) as report:
            reader = csv.reader(report, delimiter='\t')
            for line in reader:
                if line[1] == '1:115256669G>A':
                    if 'NM_001007553' in line[29]:
                        self.assertEqual(line[3], 'True')


    def test_preferred_transcripts_high_strictness_false(self):
        """
        Check that preferred transcripts are labelled correctly. 
        XM_005245221 should be false with high strictness and true with
        low strictness.
        """
        self.pt = preferred_transcripts()
        self.pt.load(os.path.abspath('test/PreferredTranscripts.txt'))
        self.pt.apply(self.report.report_path, 'high')
        with open(self.report.report_path) as report:
            reader = csv.reader(report, delimiter='\t')
            for line in reader:
                if line[1] == '1:162681151T>G':
                    if 'XM_005245221' in line[29]:
                        self.assertEqual(line[3], 'False')


    def test_preferred_transcripts_low_strictness_true_1(self):
        """
        Check that preferred transcripts are labelled correctly. 
        NM_001007553 should be true with high and low strictness
        """
        self.pt = preferred_transcripts()
        self.pt.load(os.path.abspath('test/PreferredTranscripts.txt'))
        self.pt.apply(self.report.report_path, 'low')
        with open(self.report.report_path) as report:
            reader = csv.reader(report, delimiter='\t')
            for line in reader:
                if line[1] == '1:115256669G>A':
                    if 'NM_001007553' in line[29]:
                        self.assertEqual(line[3], 'True')


    def test_preferred_transcripts_low_strictness_true_2(self):
        """
        Check that preferred transcripts are labelled correctly. 
        XM_005245221 should be false with high strictness and true with
        low strictness.
        """
        self.pt = preferred_transcripts()
        self.pt.load(os.path.abspath('test/PreferredTranscripts.txt'))
        self.pt.apply(self.report.report_path, 'low')
        with open(self.report.report_path) as report:
            reader = csv.reader(report, delimiter='\t')
            for line in reader:
                if line[1] == '1:162681151T>G':
                    if 'XM_005245221' in line[29]:
                        self.assertEqual(line[3], 'True')


    def test_preferred_transcripts_false(self):
        """
        Check that preferred transcripts are labelled correctly. 
        Not in preferred transcripts - should be false.
        """
        self.pt = preferred_transcripts()
        self.pt.load(os.path.abspath('test/PreferredTranscripts.txt'))
        self.pt.apply(self.report.report_path, 'low')
        with open(self.report.report_path) as report:
            reader = csv.reader(report, delimiter='\t')
            for line in reader:
                if line[1] == '3:41265953CT>C':
                    self.assertEqual(line[3], 'False')


    # apply bed files
    def test_bed_files(self):
        """
        Check that correct number of variants are filtered out with BED file
        Should be 6
        """
        self.bed = bed_object()
        self.bed.apply_single(
            os.path.abspath('test/test_bed_files/bed1.bed'), self.report
            )
        n = sum(1 for line in open(os.path.abspath(
            'test/test_bed_files/SAMPLE1_bed1_VariantReport.txt'
            )))
        self.assertEqual(n , 6)


    def test_bed_files_multiple(self):
        """
        Check that correct number of variants are filtered out with 
        multiple BED files
        Should be 6 and 10
        """
        self.bed = bed_object()
        self.bed.apply_multiple(
            os.path.abspath('test/test_bed_files/'), self.report
        )
        bed1_sum = sum(1 for line in open(os.path.abspath(
            'test/test_bed_files/SAMPLE1_bed1_VariantReport.txt'
            )))
        bed2_sum = sum(1 for line in open(os.path.abspath(
            'test/test_bed_files/SAMPLE1_bed2_VariantReport.txt'
            )))
        self.assertEqual(bed1_sum, 6)
        self.assertEqual(bed2_sum, 10)


# Runs all tests when the script is run
# command: python -m unittest -v test
if __name__ == '__main__':
    unittest.main()
