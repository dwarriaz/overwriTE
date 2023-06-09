import sys
import csv

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
    

def main():
    
    gene_list_reader=readCSV(sys.argv[1])
    
    gene_list=gene_list_reader.csvreader()

    batchlist = []
    
    step = len(gene_list)//36
    for x in range(0,len(gene_list),step):
        package = gene_list[x:x+step]
        batchlist.append(package)

    counter = 1
    for batch in batchlist:
        final_format='('
        for gene in batch:
            final_format += "'"+gene[0]+"'"+","
        final_format = final_format.rstrip(final_format[-1])
        final_format += ")"
        
        with open('output/Batch'+str(counter),'w') as file:
            file.write(final_format)
        counter += 1

    

if __name__ == "__main__":
    main()
    
