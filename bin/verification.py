import csv
from Bio import SeqIO


filename = 'verifiedinteractions.tsv'
geneDir = 'query'
column_name = 'mutation-sRNA&gene' # 6
fileformat = ".fasta"



rows = []
column_names = []
organismFolder = []
mutation = []
geneFile = []


def parseFiletsv():
    with open(filename) as tsvfile:
        tsvreader = csv.reader(tsvfile, delimiter="\t")
        for index, dataLine in enumerate(tsvreader):
            if(index == 0):
                column_names.append(dataLine) # For column Names.
            else:
                rows.append(dataLine)



def getMutationgene():
    totalOccurances = 0
    for index, element in enumerate(rows):
        if(len(element) > 6):
            if element[6] != '':
                if( '@' not in element[6]): # To remove data where it says @martin.
                    geneFile.append(element[0])
                    organismFolder.append(element[3])
                    mutation.append(element[6])
                    totalOccurances += 1

    print("Total "+ str(totalOccurances) + " found in " + filename + " file")





def process():
    for index, element in enumerate(geneFile):
        completeFileName = element+"_"+organismFolder[index]+fileformat
        fasta_sequences = SeqIO.parse(open("input/"+organismFolder[index]+"/query/"+completeFileName),'fasta')
        for fasta in fasta_sequences:
            fastaSequence = fasta.seq
            fastaSeqCharArray = str(fastaSequence).split(",") # Convert to string array for checkig easy purpose. for index wise search.
            print("Fasta Sequence : \n", fastaSeqCharArray)
            print("To Check mutation:", mutation[index], "\n \n \n")


parseFiletsv()
getMutationgene()
process()
