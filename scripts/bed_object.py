import os
import csv


# set global variable to point to bedtools executables
BEDTOOLS_PATH = '/Users/erik/Applications/bedtools2/bin/'


# ----------------- BED CLASS -----------------------------------------

class bed_object:
    def make_intersect_bed(self, bedfile, in_vcf):
        # turn report into bed file and save as temp.bed file
        with open(in_vcf.report_path, 'r') as report:
            reader = csv.reader(report, delimiter='\t')
            with open(in_vcf.report_path + '.temp', 'w') as out: 
                for line in reader:
                    if line[0][0] != '#':
                        variant = line[1].split(':')
                        out.write(variant[0] + "\t"
                            + str(int(variant[1].strip('AGTC>,')) - 1)
                            + "\t" + variant[1].strip('AGTC>,') + "\t"
                            + line[1] + "\n")
            out.close()
        report.close()

        # make output folder if it doesnt exist, based on input folder name
        in_folder = os.path.split(os.path.dirname(bedfile))[-1]
        out_folder = os.path.join(in_vcf.output_dir, in_folder)
        if not os.path.exists(out_folder):
            os.mkdir(out_folder)

        # intersect report bed and input bed, save as variable
        report_bed = in_vcf.report_path + '.temp'
        intersect_command = (BEDTOOLS_PATH + 'intersectBed -a ' + report_bed + ' -b ' + bedfile)
        try:
            results_intersect = os.popen(intersect_command).read()
        except IOError:
            print('BEDTools error')

        # write bed variable to file, remove report bed
        self.bed_name = os.path.basename(bedfile).rstrip('.bed')
        self.intersect_bed = os.path.join(out_folder, (in_vcf.sample + '_' + self.bed_name + '_intersect.bed'))
        self.out_folder = out_folder

        out = open(self.intersect_bed, 'w') 
        out.write(results_intersect)
        out.close()
        os.remove(report_bed)


    def apply_bed(self, bedfile, in_vcf):
        # load temp intersect bed, read variant into list
        with open(self.intersect_bed, 'r') as bed:
            results = csv.reader(bed, delimiter='\t')
            keep = []
            for line in results:
                keep.append(line[3])

        # open empty file
        outfile = os.path.join(self.out_folder, (in_vcf.sample + '_' + self.bed_name + '_VariantReport.txt'))
        bed = open(outfile, 'w')
        bed_report = csv.writer(bed, delimiter='\t')

        # loops through original report, if there is a match with the bed list, keep the line
        with open(in_vcf.report_path) as report:
            results = csv.reader(report, delimiter='\t')
            for line in results:
                if line[0][0] == '#':
                    bed_report.writerow(line)
                if line[1] in keep:
                    bed_report.writerow(line)
        
        print('Out: ' + outfile)
        os.remove(os.path.abspath(self.intersect_bed))


    def apply_single(self, bedfile, in_vcf):
        self.make_intersect_bed(bedfile, in_vcf)
        self.apply_bed(bedfile, in_vcf)

    
    def apply_multiple(self, bed_folder, in_vcf):
        for file in os.listdir(bed_folder):
            in_bed = os.path.join(bed_folder, file)
            print('In:  ' + in_bed)
            self.make_intersect_bed(in_bed, in_vcf)
            self.apply_bed(in_bed, in_vcf)
