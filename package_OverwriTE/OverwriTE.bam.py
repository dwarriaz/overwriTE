#! usr/bin/python3
import sys

def transcriptBuilder(samfile):
    pairedDict = {} 
    with open(samfile) as sFile:
        for line in sFile:
            data = line.strip('\n').split('\t')
            if data[0] in pairedDict.keys():
                stop = int(data[3]) - int(data[8])
                pairedDict[data[0]].append(stop)
            else:
                pairedDict[data[0]] = [data[2],int(data[3])]
    return pairedDict

def main():

   # for entry in transcriptBuilder(sys.argv[1]).values():
   print(transcriptBuilder(sys.argv[1]))


if __name__ == "__main__":
    main()