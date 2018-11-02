# vcf_parser
Code for parsing VCF files and making a text report.

Special formatting:  
- Allele frequency has been calculated from the AD
- Genotypes are reformatted from 0/0, 0/1 and 1/1 to HOM_REF, HET and HOM_ALT, respectively, to prevent them from appearing as dates in Excel
- dbSNP, COSMIC and HGMD are parsed from the 'Existing_variation' VEP field
- All population frequencies from ExAC and 1KG are re-formatted to a percentage
- Intron and exon numberings are renamed from x/y to x|y (where x and y are numbers), to prevent them appearing as dates in Excel
- HGVS coding (HGVSc) and protein (HGVSp) sequences have had the transcript name trimmed off



