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
                    # print("Negative Position Encountered in sRNA check")
                    print()
                elif("-" in geneCheck):
                    # print("Negative Position encountered in genecheck")
                    print()
                else:
                    boolsRNA(sRNACheck, sRNAGene[record].get("fasta") ,sRNAGene[record].get("name"))
                    boolsGene(geneCheck, sRNAGene[record].get("gene").get("Sequence") ,sRNAGene[record].get("gene").get("name"))
                    complementBase(sRNACheck, geneCheck)

        elif ("," in mutation):
            print("\n Mutation with Bar ',' : \t", mutation)
        else:
            geneCheck, sRNACheck = mutation.split('&', 1)
            if("-" in sRNACheck):
                # print("Negative Position Encountered in sRNA check")
                print()
            elif("-" in geneCheck):
                # print("Negative Position encountered in genecheck")
                print()
            else:
                boolsRNA(sRNACheck, sRNAGene[record].get("fasta") ,sRNAGene[record].get("name"))
                boolsGene(geneCheck, sRNAGene[record].get("gene").get("Sequence") ,sRNAGene[record].get("gene").get("name"))
                complementBase(sRNACheck, geneCheck)


def boolsRNA(checkMut, sequence, name):
    position = int(re.findall(r'\d+', checkMut)[0])
    base_1 = checkMut[0].upper()
    base_2 = checkMut[-1].upper()
    if base_1 == "U":
        base_1 = "T"
    if len(sequence) > position:
        fromSeq = sequence[position-1].upper()
        if fromSeq == base_1:
            print()
            # print("\n  To check ",base_1," in position ",position," in from sequence in file ",name," : Evaluation result Got ",fromSeq," : True")
            # complementBase(base_1, base_2)
        else :
            print("\n  To check ",base_1," in position ",position," in from sequence in file ",name," : Evaluation result Got ",fromSeq," : False")
            # complementBase(base_1, base_2)
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
            print()
            # print("\n  To check ",base_1," in position ",position," in from sequence in file ",name," : Evaluation result Got ",fromSeq," : True")
            # complementBase(base_1, base_2)
        else :
            print("\n  To check ",base_1," in position ",position," in from sequence in file ",name," : Evaluation result Got ",fromSeq," : False")
            # complementBase(base_1, base_2)
    else :
        print("Check for Typo as Sequence is shorter than specified index to check in file", name)
        print("Lenght of Sequence : \t",len(sequence))
        print("Encountered position to find : \t", position)

""" Watson-crick Method for finding Gene complement """
def complementBase(before, after):
    print(before, after)
    b_check_1 = before[0].upper()
    a_check_1 = after[0].upper()
    b_check_2 = before[-1].upper()
    a_check_2 = after[-1].upper()
    check = [b_check_1, a_check_1, b_check_2 , a_check_2]
    G = "G"
    C = "C"
    A = "A"
    U = "U"
    T = "T"

    for element in check:
        if element == T:
            element = U

    if ((b_check_1 == C and a_check_1 == G) or (b_check_1 == G and a_check_1 == C)):
        print("Gene complement matched", b_check_1, " <-> ", a_check_1)
    elif ((b_check_1 == A and a_check_1 == U) or (b_check_1 == U and a_check_1 == A)):
        print("Gene complement matched", b_check_1, " <-> ", a_check_1)
    else:
        print("Gene complement mismatched", b_check_1, " <-!WRONG!-> ", a_check_1)

    if ((b_check_2 == C and a_check_2 == G) or (b_check_2 == G and a_check_2 == C)):
        print("Gene complement matched", b_check_2, " <-> ", a_check_2)
    elif ((b_check_2 == A and a_check_2 == U) or (b_check_2 == U and a_check_2 == A)):
        print("Gene complement matched", b_check_2, " <-> ", a_check_2)
    else:
        print("Gene complement mismatched", b_check_2, " <-!WRONG!-> ", a_check_2)


""" Run this only when you want to update the data  """

parseFiletsv()
getMutationgene()
getSRNA()
getGene()
saveFile()



""" Read existing data from json file. """
# readData()

""" Entry point for the Application. """
# processLoop()


