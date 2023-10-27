import sys
import csv
import pandas as pd 
import concurrent.futures


class CommandLine:
    
    def __init__(self):
        
        import argparse
        
        self.parser = argparse.ArgumentParser()
        
        self.parser.add_argument('-Gene', help = 'Path to csv containing ranges of Exons and Introns')
        self.parser.add_argument('-Output', help = 'Path to csv for overwriTE results')
        
        self.args = self.parser.parse_args()

class readCSV:
    def __init__(self, filename= None):
        self.filename = filename
    
    def csvreader(self):
        csvdata = [] 
        with open(self.filename,'r') as file:
            csvreader = csv.reader(file)
            for row in csvreader:
                csvdata.append(row)
        return csvdata

class table_parser:
    def __init__(self):
        self.chrom_library= {'chr1':[],'chr2':[],'chr3':[],'chr4':[],'chr5':[],'chr6':[],
        'chr7':[],'chr8':[],'chr9':[],'chr10':[],'chr11':[],'chr12':[],'chr13':[],'chr14':[],'chr15':[],'chr16':[],
        'chr17':[],'chr18':[],'chr19':[],'chr20':[],'chr21':[],'chr22':[],'chrX':[],'chrY':[]}
    
    def table_split(self,rmsk_table):
        for entry in rmsk_table:
            if entry[0] in self.chrom_library.keys():
                self.chrom_library[entry[0]].append(entry)
        return self.chrom_library
    
    def table_sort(self):
        for repeats in self.chrom_library.values():
            repeats.sort(key = lambda repeats: repeats[1])
        return self.chrom_library

class OverwriTE_Utility:
    
    def __init__(self,sorted_chrom_library):
       self.chrom_library = sorted_chrom_library
        
    def binary_search(self,chrom,rmsk_index):
        for repeats in self.chrom_library[chrom]:

            low = 0
            high = (len(self.chrom_library[chrom])-1)

            while low <= high:
                mid = (high + low)//2
                
                if int(rmsk_index) > int(self.chrom_library[chrom][mid][1]):
                    low = mid+1
                    if rmsk_index in range(int(self.chrom_library[chrom][mid][1]),int(self.chrom_library[chrom][mid][2])):
                        return mid
                elif int(rmsk_index) < int(self.chrom_library[chrom][mid][1]):
                    high = mid - 1
                    if int(rmsk_index) in range(int(self.chrom_library[chrom][mid][1]),int(self.chrom_library[chrom][mid][2])):
                        return mid
                else: 
                    return mid
    
    def table_pull(self, chrom, index1, index2):
        #for repeats in self.chrom_library[chrom]:
        package = [self.chrom_library[chrom][index1:index2]]
        return package 
    
    def trim(arg):
        arg.strip('"')
        temp_list = arg.split(",")
        temp_list.pop(-1)
        return temp_list
    
    def Exon_rangelist_gen(uniq_start_list, uniq_end_list):
        rangelist = []
        for x in range(0,len(uniq_start_list)):
            package = range(int(uniq_start_list[x]),int(uniq_end_list[x]))
            rangelist.append(package)
        return rangelist
    

class Isoform:
    def __init__(self, isoform):
        self.name = isoform[0]
        self.name2 = isoform[1]
        if len(isoform) == 
        self.chrom = if len(isoform[2]) 
        self.strand = isoform[3]
        self.txStart = isoform[4]
        self.txEnd = isoform[5]
        self.cdsStart = isoform[6]
        self.cdsEnd = isoform[7]
        self.exonCount = isoform[8]
        self.exonStarts = isoform[9]
        self.exonEnds = isoform[10]

        self.repeats = list()

    def rmsk_finder(self):
        
        index_for_start = int()
        index_for_stop = int()
        
        if self.strand == '+':
            index_for_start = util.binary_search(self.chrom,int(self.txStart)-1500)
            index_for_stop = util.binary_search(self.chrom,int(self.txEnd)+300)
        elif self.strand == '-':
            index_for_start = util.binary_search(self.chrom,int(self.txStart)-300)
            index_for_stop = util.binary_search(self.chrom,int(self.txEnd)+1500)
        
        self.repeats = util.table_pull(self.chrom,index_for_start,index_for_stop)
    
    def txStart_UTR(self):
        output_list = [] 
        if self.strand == '+':
            for repeats in self.repeats:
                for repeat in repeats:
                    repeat_set = set(range(int(repeat[1]),int(repeat[2])))
                
                    if len(repeat_set.intersection(set(range((int(self.txStart) - 1500),int(self.txStart))))) !=0:
                        package = [self.name,self.name2,self.chrom,self.strand,int(self.txStart),
                            int(self.txEnd), repeat[4], repeat[5], repeat[6], repeat[7],'Promoter_Region', 
                            len(repeat_set), len(repeat_set.intersection(set(range((int(self.txStart) - 1500),int(self.txStart)))))]

                        output_list.append(package)
                
                    if len(repeat_set.intersection(set(range(int(self.txEnd),int(self.txEnd) + 300)))) !=0:
                        package = [self.name,self.name2,self.chrom,self.strand,int(self.txStart),
                            int(self.txEnd), repeat[4], repeat[5], repeat[6], repeat[7],'3_UTR', 
                            len(repeat_set), len(repeat_set.intersection(set(range(int(self.txEnd),int(self.txEnd) + 300))))]

                        output_list.append(package)
        
        return output_list
                

        
    def Overlap(self):
        start_list = OverwriTE_Utility.trim(self.exonStarts)
        end_list = OverwriTE_Utility.trim(self.exonEnds)

        exon_range_list = OverwriTE_Utility.Exon_rangelist_gen(start_list,end_list)

        output_list = []

        for exon in exon_range_list:
            for repeats in self.repeats:
                for repeat in repeats: 
                    exon_set = set(exon)

                    repeat_set = set(range(int(repeat[1]),int(repeat[2])))

                    if len(exon_set.intersection(repeat_set)) >= 1:
                        if len(exon_set) == len(repeat_set):
                            package = [self.name,self.name2,self.chrom,self.strand,int(self.txStart),
                            int(self.txEnd), repeat[4], repeat[5], repeat[6], repeat[7],'exon', 
                            len(repeat_set), len(exon_set.intersection(repeat_set))]
                            output_list.append(package)
                        elif len(exon_set) != len(repeat_set):
                            package = [self.name,self.name2,self.chrom,self.strand,int(self.txStart),
                            int(self.txEnd), repeat[4], repeat[5], repeat[6], repeat[7],'junction', 
                            len(repeat_set), len(exon_set.intersection(repeat_set))]
                            output_list.append(package)
                    elif len(exon_set.intersection(repeat_set)) == 0:
                        package = [self.name,self.name2,self.chrom,self.strand,int(self.txStart),
                        int(self.txEnd), repeat[4], repeat[5], repeat[6], repeat[7],'intron', 
                        len(repeat_set), len(exon_set.intersection(repeat_set))]
                        output_list.append(package)
        
        return output_list

class Parallel_overwriTE(): 
    def execute(isoform):
        work = Isoform(isoform)
        work.rmsk_finder()
        output = work.Overlap()
        second_output = work.txStart_UTR()
        output_package = output + second_output

        return output_package



    


def main(util,transcript_table):
    
    #ThisCommandLine = CommandLine()

    ##Creates CSV reading object and executes 

    '''rmsk_reader = readCSV('output/rmsk_table.csv')
    rmsk_table = rmsk_reader.csvreader()
    
    ## Will split the rmsk table by chromosome 

    parser=table_parser()
    chrom_library=parser.table_split(rmsk_table)

    ## Will return a dictionary whose key is the chrom and contains a lists of lists 
    
    sorted_chrom_library=parser.table_sort()'''

    ## Read in CSV for Transcription Data 

    '''transcript_reader = readCSV(ThisCommandLine.args.Gene) 
    transcript_table = transcript_reader.csvreader()'''

    ## Begining OverwriTE functionality
    
    output = []
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
        results = executor.map(Parallel_overwriTE.execute, transcript_table)
        for result in results:
            output += result 

    overlaps = pd.DataFrame(output)
    overlaps.to_csv(ThisCommandLine.args.Output, index=False, mode='a')
    

        


if __name__ == "__main__":
    ThisCommandLine = CommandLine()
    rmsk_reader = readCSV('output/rmsk_table.csv')
    rmsk_table = rmsk_reader.csvreader()
    parser=table_parser()
    chrom_library=parser.table_split(rmsk_table)
    sorted_chrom_library=parser.table_sort()
    util = OverwriTE_Utility(sorted_chrom_library)

    transcript_reader = readCSV(ThisCommandLine.args.Gene) 
    transcript_table = transcript_reader.csvreader()
    print('rmsk_sort completed')
    
    main(util,transcript_table)

        