import os
import csv
import warnings


# ----------------- PREFERRED TRANSCRIPTS CLASS -----------------------

class preferred_transcripts:
    def load(self, args):
        # Load in external file containing preferred transcripts
        try:
            preferred_transcripts = []
            with open(args.transcripts, 'r') as pt:
                pt_reader = csv.reader(pt, delimiter='\t')
                for line in pt_reader:
                    preferred_transcripts += [line[1]]
            self.list = preferred_transcripts
        except:
            warnings.warn('Could not find preferred transcripts file, continuing without.')


    def apply(self, report, transcript_id='Feature'):
        # open report file and empty temp file to save new output
        report_temp = os.path.join(report + '.temp')
        f1 = open(report, 'rb')
        reader = csv.reader(f1, delimiter='\t')
        f2 = open(report_temp, 'wb')
        writer = csv.writer(f2, delimiter='\t')

        # add header to new file, find the column that the transcript ID is kept in
        header = next(reader)
        writer.writerow(header)
        preferred_column = None
        transcript_column = None
        preferred_column = header.index('Preferred')
        transcript_column = header.index(transcript_id)

        # if transcript ID cant be found, exit the script and tidy up. The rest of the functions can carry on as normal becuase the originial variant report hasnt been touched.
        if transcript_column == None or preferred_column == None:
            warnings.warn('Could not find transcripts/ preferred column in variant report file, continuing without adding preferred transcripts.')
            f1.close()
            f2.close()
            os.remove(report_temp)
            return

        # loop though the report, change preferred to true if there is a match, otherwise just copy the row
        for row in reader:
            try:
                if row[transcript_column-1] in self.list:
                    writer.writerow(row[0:preferred_column-1] + ['True'] + row[preferred_column:])
                else:
                    writer.writerow(row)
            except:
                writer.writerow(row)

        # tidy up
        f1.close()
        f2.close()
        os.remove(report)
        os.rename(report_temp, report)