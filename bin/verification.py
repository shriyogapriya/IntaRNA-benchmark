import csv
from Bio import SeqIO
import json
import re


filename = 'verified_interactions.tsv'
column_name = 'mutation-sRNA&gene' # 6
fileformat = ".fasta"



rows = []
column_names = []
organismFolder = []
mutation = []
geneFile = []
sRNAFile = []


sRNAGene = {}

""" parses the tsv file and gets the tabular data inside an arraz named <rows> """
def parseFiletsv():
    with open(filename) as tsvfile:
        tsvreader = csv.reader(tsvfile, delimiter="\t")
        for index, dataLine in enumerate(tsvreader):
            if(index == 0):
                column_names.append(dataLine) # For column Names.
            else:
                rows.append(dataLine)

""" Segregate the data obtained from the tsv file"""
def getMutationgene():
    totalOccurances = 0
    for index, element in enumerate(rows):
        if(len(element) > 6):
            if element[6] != '':
                if( '@' not in element[6]): 
                    sRNAFile.append(element[0])
                    geneFile.append(element[1])
                    organismFolder.append(element[3])
                    mutation.append(element[6])
                    totalOccurances += 1
    print("Total "+ str(totalOccurances) + " found in " + filename + " file")



""" gets the mutation gene file in the given geneDirectory """
def getGene():
    for index, element in enumerate(geneFile):
        fasta_sequences = SeqIO.parse(open("input/"+organismFolder[index]+"/target/"+organismFolder[index]+".fa"),'fasta')
        for Jindex, fasta in enumerate(fasta_sequences):
            if(element == fasta.name):
                fastaSequence = fasta.seq
                fastaSeqCharArray = str(fastaSequence).split(",")[0] # Convert to string array for checkig easy purpose. for index wise search.
                sRNAGene[index]["gene"]= {"name": fasta.name , "Sequence": fastaSeqCharArray}
                # print("Gene Name  : \n", fasta.name)
                # print("Fasta Sequence : \n", fastaSeqCharArray)



""" gets the mutation from sRNA file in the given sRNADirectory """
def getSRNA():
    for index, element in enumerate(sRNAFile):
        completeFileName = element+"_"+organismFolder[index]+fileformat
        fasta_sequences = SeqIO.parse(open("input/"+organismFolder[index]+"/query/"+completeFileName),'fasta')
        for fasta in fasta_sequences:
            fastaSequence = fasta.seq
            fastaSeqCharArray = str(fastaSequence).split(",")[0] # Convert to string array for checkig easy purpose. for index wise search.
            sRNAGene[index] = {"fasta": fastaSeqCharArray,  "mutation":mutation[index], "name":completeFileName }



""" Save final sRNAGENE data file """
def saveFile():
    with open('sRNAGeneComplete.json', 'w') as outfile:  
        json.dump(sRNAGene, outfile)

""" Get Data from file fropm previous run """
def readData():
    with open('sRNAGeneComplete.json') as json_file:  
        retThis = json.load(json_file)
        print("Data Imported Succesfully!")
        return retThis

    


""" Main process loop for checking and evaluating based on mutation"""
def processLoop():
    sRNAGene = readData()
    for index, record in enumerate(sRNAGene):
        mutation = sRNAGene[record].get("mutation")
        if ("|" in mutation):
            print("\n Multiple Evaluations using : \t", mutation)
            multiEval = mutation.split("|")
            for index , element in enumerate(multiEval):
                geneCheck, sRNACheck = element.split('&', 1)
                if("-" in sRNACheck):
                    print("Negative Position Encountered in sRNA check")
                elif("-" in geneCheck):
                    print("Negative Position encountered in genecheck")
                else:
                    boolsRNA(sRNACheck, sRNAGene[record].get("fasta") ,sRNAGene[record].get("name"))
                    boolsGene(geneCheck, sRNAGene[record].get("gene").get("Sequence") ,sRNAGene[record].get("gene").get("name"))
        elif ("," in mutation):
            print("\n Mutation with Bar ',' : \t", mutation)
        else:
            geneCheck, sRNACheck = mutation.split('&', 1)
            if("-" in sRNACheck):
                print("Negative Position Encountered in sRNA check")
            elif("-" in geneCheck):
                print("Negative Position encountered in genecheck")
            else:
                boolsRNA(sRNACheck, sRNAGene[record].get("fasta") ,sRNAGene[record].get("name"))
                boolsGene(geneCheck, sRNAGene[record].get("gene").get("Sequence") ,sRNAGene[record].get("gene").get("name"))


def boolsRNA(checkMut, sequence, name):
    position = int(re.findall(r'\d+', checkMut)[0])
    base_1 = checkMut[0].upper()
    base_2 = checkMut[-1].upper()
    if base_1 == "U":
        base_1 = "T"
    if len(sequence) > position:
        fromSeq = sequence[position-1].upper()
        if fromSeq == base_1:
            print("\n  To check ",base_1," in position ",position," in from sequence in file ",name," : Evaluation result Got ",fromSeq," : True")
            complementBase(base_1, base_2)
        else :
            print("\n  To check ",base_1," in position ",position," in from sequence in file ",name," : Evaluation result Got ",fromSeq," : False")
            complementBase(base_1, base_2)
    else :
        print("Check for Typo as Sequence is shorter than specified index to check in file", name)
        print("Lenght of Sequence : \t",len(sequence))
        print("Encountered position to find : \t", position)

def boolsGene(checkMut, sequence, name):
    position = int(re.findall(r'\d+', checkMut)[0])
    base_1 = checkMut[0].upper()
    base_2 = checkMut[-1].upper()
    if base_1 == "U":
        base_1 = "T"
    if len(sequence) > position:
        fromSeq = sequence[position-1].upper()
        if fromSeq == base_1:
            print("\n  To check ",base_1," in position ",position," in from sequence in file ",name," : Evaluation result Got ",fromSeq," : True")
            complementBase(base_1, base_2)
        else :
            print("\n  To check ",base_1," in position ",position," in from sequence in file ",name," : Evaluation result Got ",fromSeq," : False")
            complementBase(base_1, base_2)
    else :
        print("Check for Typo as Sequence is shorter than specified index to check in file", name)
        print("Lenght of Sequence : \t",len(sequence))
        print("Encountered position to find : \t", position)

""" Watson-crick Method for finding Gene complement """
def complementBase(base_1, base_2):
    G = "G"
    C = "C"
    A = "A"
    U = "U"
    T = "T"

    if base_1 == T:
        base_1 =  U
    elif base_2 == T:
        base_2 =  U

    if ((base_1 == C and base_2 == G) or (base_1 == G and base_2 == C)):
        print("Gene complement matched", base_1, " <-> ", base_2)
    elif ((base_1 == A and base_2 == U) or (base_1 == U and base_2 == A)):
        print("Gene complement matched", base_1, " <-> ", base_2)
    else:
        print("Gene complement mismatched", base_1, " <-> ", base_2)


""" Run this only when you want to update the data  """

# parseFiletsv()
# getMutationgene()
# getSRNA()
# getGene()
# saveFile()



""" Read existing data from json file. """
# readData()

""" Entry point for the Application. """
processLoop()


