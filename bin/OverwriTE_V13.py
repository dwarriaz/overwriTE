import csv
import pandas as pd 
import argparse

class CommandLine:
    
    def __init__(self):
        
        self.parser = argparse.ArgumentParser()
        
        self.parser.add_argument('-SQL', help = 'Path to SQL query')
        self.parser.add_argument('-Output', help = 'File name for output')
        
        self.args = self.parser.parse_args()

class readCSV:
    
    def __init__(self, filename= None):
        self.filename = filename
    
    def csvreader(self):
        csvdata = [] 
        with open(self.filename,'r') as file:
            csvreader = csv.reader(file,delimiter="\t")
            for row in csvreader:
                csvdata.append(row)
        return csvdata

class OverwriTE:
    
    def __init__(self,row):
        
        self.name = row[0]
        self.name2 = row[1]
        self.chrom = row[2]
        self.genoStrand = row[3]
        self.txStart = int(row[4])
        self.txEnd = int(row[5])
        self.exonCount = row[6]
        self.exonStarts = (row[7])
        self.exonEnds = (row[8])
        self.repStart = int(row[9])
        self.repEnd = int(row[10])
        self.repName = row[11]
        self.repFamily = row[12]
        self.repClass = row[13]
        self.repStrand = row[14]
        self.cdsStart = int( row[15])
        self.cdsEnd = int(row[16])

    def package(self,classification,len_classification,region_length):
        complete_seq_name = self.repName +'_range='+self.chrom+':'+str(self.repStart)+'-'+str(self.repEnd)+'_strand='+self.repStrand

        OverwriTE_entry = [self.name, self.name2, self.chrom, self.genoStrand,(self.txEnd-self.txStart), self.repName, 
                           self.repStart,(self.repEnd-self.repStart),self.repFamily, self.repClass, self.repStrand, classification,len_classification,
                           region_length,complete_seq_name]
        #print(complete_seq_name,classification,region_length)
        return OverwriTE_entry

    def classification(self): 

        repSet = set(range(self.repStart,self.repEnd,1))

        if self.genoStrand == '+':
            if self.repEnd in range(self.txStart-1500,self.txStart):
                return self.package('Promoter_Region',len(repSet),1500)
            elif self.repStart in range(self.txEnd,(self.txEnd + 300),1):
                return self.package ('3_End_Region',len(repSet),300)
        elif self.genoStrand == '-':
            if self.repStart in range(self.txEnd,self.txEnd+1500):
                return self.package('Promoter_Region',len(repSet),1500)
            elif self.repEnd in range((self.txStart-300),self.txStart,1):
                return self.package('3_End_Region',len(repSet),300)

        exonStartlist = self.exonStarts.split(',')
        exonStoplist = self.exonEnds.split(',')

        exonStartlist.pop(-1)
        exonStoplist.pop(-1)
        
        totalExonLength = int()
        for x in range(0,len(exonStartlist)):
            totalExonLength += (int(exonStoplist[x])-int(exonStartlist[x]))

        for x in range(0,len(exonStartlist)):
            calculation = repSet.intersection(range(int(exonStartlist[x]),int(exonStoplist[x])))

            if len(calculation) == 0:
                genoRegion = (self.txEnd - self.txStart)-totalExonLength
                return self.package('Intron',len(repSet),genoRegion)

            elif len(calculation) != len(repSet):
                if self.genoStrand == "+":
                    if self.txStart in repSet:
                        return self.package('Transcription_Start_Site',self.txStart,1) 
                    elif self.txEnd in repSet:
                        return self.package('Transcription_End_Site',self.txEnd,1)
                elif self.genoStrand == "-":
                    if self.txStart in repSet:
                        return self.package('Transcription_End_Site',self.txStart,1) 
                    elif self.txEnd in repSet:
                        return self.package('Transcription_Start_Site',self.txEnd,1)

                if int(exonStartlist[x]) in repSet: 
                    junctionPosition = exonStartlist[x]
                    return self.package('Junction',junctionPosition, len(calculation))
                elif int(exonStoplist[x]) in repSet:
                    junctionPosition = exonStoplist[x]
                    return self.package('Junction',junctionPosition, len(calculation))
            
            elif len(calculation) == len(repSet):
                if self.txStart == self.cdsStart and self.cdsEnd == self.cdsStart:
                    return self.package('Exon',len(repSet),totalExonLength)
                
                else:
                    if self.genoStrand == '+':
                        if self.repStart <= self.cdsStart:
                            return self.package('5UTR', len(calculation),(self.cdsStart-self.txStart))
                        if self.repStart >= self.cdsEnd:
                            return self.package('3UTR', len(calculation),(self.txEnd-self.cdsEnd))
                    elif self.genoStrand == '-':
                        if self.repStart <= self.cdsStart:
                            return self.package('3UTR', len(calculation),(self.cdsStart-self.txStart))
                        if self.repStart >= self.cdsEnd:
                            return self.package('5UTR', len(calculation),(self.txEnd-self.cdsEnd))
                    return self.package('Exon',len(repSet),(totalExonLength))



def main():
    
    ThisCommandLine = CommandLine()

    print('Creating:'+ThisCommandLine.args.Output)
    
    table_reader = readCSV(ThisCommandLine.args.SQL)
    table = table_reader.csvreader()
    table.pop(0)

    result = [OverwriTE(row).classification() for row in table]
    
    overlaps = pd.DataFrame(result, columns = ['name', 'name2', 'chrom', 'genoStrand','genoLength', 'repName','repStart','repLength', 'repFamily', 'repClass', 'repStrand', 'classification', 'len_classification','regionLength','compleTE_seq'])
    overlaps.to_csv(ThisCommandLine.args.Output, index = False)
    

if __name__ == "__main__":
    main()