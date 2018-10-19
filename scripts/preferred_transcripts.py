import os
import csv


# -- PREFERRED TRANSCRIPTS CLASS --------------------------------------

class preferred_transcripts:
    
    def load(self, args):
        '''
        Take a preferred transcript file and extract the preferred 
        transcripts into a list, and save as self.list. This function is 
        only called if a preferred transcripts file is included, so a 
        warning is thrown if it can't find the file.
        '''
        try:
            preferred_transcripts = []
            with open(args.transcripts, 'r') as pt:
                pt_reader = csv.reader(pt, delimiter='\t')
                for line in pt_reader:
                    preferred_transcripts += [line[1]]
            self.list = preferred_transcripts
        except:
            self.list = None


    def apply(self, report, transcript_id='Feature'):
        '''
        Take a variant report and loop through each row, change 
        preferred transcript to true if there's a match, otherwise 
        change to false.
        '''
        if self.list:
            # open report file and empty temp file to save new output
            report_temp = os.path.join(report + '.temp')
            f1 = open(report, 'rb')
            reader = csv.reader(f1, delimiter='\t')
            f2 = open(report_temp, 'wb')
            writer = csv.writer(f2, delimiter='\t')

            # add header to new file, find the column that the transcript 
            # ID is kept in
            header = next(reader)
            writer.writerow(header)
            preferred_column = None
            transcript_column = None
            preferred_column = header.index('Preferred')
            transcript_column = header.index(transcript_id)

            # if transcript ID cant be found, exit the script and tidy up.
            # The rest of the functions can carry on as normal becuase the 
            # originial variant report hasnt been touched.
            if transcript_column == None or preferred_column == None:
                print(
                    '''WARN\tCould not find transcripts/ preferred column in variant 
                    report file, continuing without adding preferred transcripts.'''
                )
                f1.close()
                f2.close()
                os.remove(report_temp)
                return

            # loop though the report, change preferred to true if there is 
            # a match, otherwise just copy the row
            for row in reader:
                if row[transcript_column] in self.list:
                    writer.writerow(row[0:preferred_column] + ['True'] + 
                    row[preferred_column+1:])
                else:
                    writer.writerow(row[0:preferred_column] + ['False'] + 
                    row[preferred_column+1:])

            # tidy up
            print('INFO\tAdded preferred transcripts.')
            f1.close()
            f2.close()
            os.remove(report)
            os.rename(report_temp, report)
        else:
            print('WARN\tCould not load preferred transcripts file provided, skipping step.')
            return