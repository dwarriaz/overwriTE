#! python3 TE_Overlap using MYSQL dual query 

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

    def package(self,classification,len_classification):
        complete_seq_name = self.repName +'_range='+self.chrom+':'+str(self.repStart)+'-'+str(self.repEnd)+'_strand='+self.repStrand

        OverwriTE_entry = [self.name, self.name2, self.chrom, self.genoStrand,(self.txEnd-self.txStart), self.repName, 
                           self.repStart,(self.repEnd-self.repStart),self.repFamily, self.repClass, self.repStrand, classification,len_classification,
                           complete_seq_name]
        return OverwriTE_entry
        
    def classification(self):
        
        repSet = set(range(self.repStart,self.repEnd,1))
        
        if self.genoStrand == '+':
            if self.repStart < self.txStart:
                return self.package('Promoter_Region',len(range(self.repStart, self.repEnd,1)))   
            
            elif self.txStart in repSet:
                return self.package('Transcription_Start_Site',self.txStart)
                
            elif self.txEnd in repSet:
                return self.package('Transcription_End_Site',self.txEnd)
            
            elif self.repStart in range(self.txEnd,(self.txEnd + 300),1):
                return self.package ('3_End_Region',len(range(self.repStart, self.repEnd,1)))

            elif self.txStart == self.cdsStart and self.cdsEnd == self.cdsStart:
                exonstartlist = self.exonStarts.split(',')
                exonstartlist.pop(-1)
                exonstoplist = self.exonEnds.split(',')
                xonstoplist.pop(-1)

                for x in range (0,len(exonstartlist)):
            
                    calculation = set(range(self.repStart,self.repEnd)).intersection(range(int(exonstartlist[x]),int(exonstoplist[x])))
                    repLen = len(range(self.repStart, self.repEnd))

                    if len(calculation) == 0:
                        return self.package('Intron', repLen)
                    elif len(calculation) > 1 and len(calculation) == repLen:
                        return self.package('Exon', len(calculation))
                    elif len(calculation) > 1 and len(calculation) != repLen: 
                        return self.package('Junction', len(calculation))
            
            else:

                exonstartlist = self.exonStarts.split(',')
                exonstartlist.pop(-1)
                exonstoplist = self.exonEnds.split(',')
                exonstoplist.pop(-1)
        

                for x in range (0,len(exonstartlist)):
            
                    calculation = set(range(self.repStart,self.repEnd)).intersection(range(int(exonstartlist[x]),int(exonstoplist[x])))
                    repLen = len(range(self.repStart, self.repEnd))

        
                    if len(calculation) == 0:
                        return self.package('Intron', repLen)
                    elif len(calculation) >1 and self.repStart < self.cdsStart:
                        return self.package('5UTR', len(calculation))
                    elif len(calculation) >1 and self.repEnd > self.cdsEnd:
                        return self.package('3UTR', len(calculation))
                    elif len(calculation) > 1 and len(calculation) == repLen:
                        return self.package('Exon', len(calculation))
                    elif len(calculation) > 1 and len(calculation) != repLen: 
                        return self.package('Junction', len(calculation))


        elif self.genoStrand == '-':

            if self.repStart > self.txEnd:
                return self.package('Promoter_Region',len(range(self.repStart, self.repEnd,1)))

            elif self.txEnd in range(self.repStart, self.repEnd,1):
                return self.package('Transcription_Start_Site',self.txEnd)

            elif self.txStart in range(self.repStart,self.repEnd):
                return self.package('Transcription_End_Site',self.txStart)

            elif self.repStart in range((self.txStart-300),self.txStart,1):
                return self.package('3_End_Region',len(range(self.repStart, self.repEnd,1)))

            elif self.txStart == self.cdsStart and self.cdsEnd == self.cdsStart:
                exonstartlist = self.exonStarts.split(',')
                exonstartlist.pop(-1)
                exonstoplist = self.exonEnds.split(',')
                exonstoplist.pop(-1)

                for x in range (0,len(exonstartlist)):
            
                    calculation = set(range(self.repStart,self.repEnd)).intersection(range(int(exonstartlist[x]),int(exonstoplist[x])))
                    repLen = len(range(self.repStart, self.repEnd))

                    if len(calculation) == 0:
                        return self.package('Intron', repLen)
                    elif len(calculation) > 1 and len(calculation) == repLen:
                        return self.package('Exon', len(calculation))
                    elif len(calculation) > 1 and len(calculation) != repLen: 
                        return self.package('Junction', len(calculation))



            else:
                exonstartlist = self.exonStarts.split(',')
                exonstartlist.pop(-1)
                exonstoplist = self.exonEnds.split(',')
                exonstoplist.pop(-1)


                for x in range (0,len(exonstartlist)):

                    calculation = set(range(self.repStart,self.repEnd)).intersection(range(int(exonstartlist[x]),int(exonstoplist[x])))
                    repLen = len(range(self.repStart, self.repEnd))


                    if len(calculation) == 0:
                        return self.package('Intron', repLen)
                    elif len(calculation) > 1 and self.repEnd > self.cdsEnd:
                        return self.package('5UTR',len(calculation))
                    elif len(calculation) > 1 and self.repStart < self.cdsStart:
                        return self.package('3UTR',len(calculation))
                    elif len(calculation) > 1 and len(calculation) == repLen:
                        return self.package('Exon', len(calculation))
                    elif len(calculation) > 1 and len(calculation) != repLen: 
                        return self.package('Junction', len(calculation))



def main():
    
    ThisCommandLine = CommandLine()

    print('Creating:'+ThisCommandLine.args.Output)
    
    table_reader = readCSV(ThisCommandLine.args.SQL)
    table = table_reader.csvreader()
    table.pop(0)
    
    result = [OverwriTE(row).classification() for row in table]
    
    overlaps = pd.DataFrame(result, columns = ['name', 'name2', 'chrom', 'genoStrand','genoLength', 'repName', 
                           'repStart','repLength', 'repFamily', 'repClass', 'repStrand', 'classification', 'len_classification','compleTE_seq'])
    overlaps.to_csv(ThisCommandLine.args.Output, index = False)
    

if __name__ == "__main__":
    main()
    
