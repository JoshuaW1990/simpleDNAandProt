import Bio
from Bio.PDB.PDBParser import PDBParser
import csv
import urllib2
import os

"""
collect the pdb id

pdbidList = []
with open("smallData.csv", 'r') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        pdbid = row[0]
        if len(pdbid) > 4:
            pdbid = pdbid[:4]
        pdbidList.append(pdbid)
pdbidList = pdbidList[1:]
result = ' '.join(pdbidList)
"""


"""
Set up a dictionary for converting the pdb 3-letters residue name to fasta 1-letter residue name
"""
pdb2fasta = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

"""
Download pdb file

print "Downlad PDB file"
pathname = "./PDB/"
count = 1
originalurl = "http://npidb.belozersky.msu.ru/data/pdb_new/biounits/1qna.pdb1.pdb"
for pdbname in pdbidList:
    print "count: ", count, "pdbname: ", pdbname
    count += 1
    filename = pdbname + ".pdb"
    completename = pathname + filename
    targeturl = originalurl.replace("1qna", pdbname)
    try:
        response = urllib2.urlopen(targeturl)
    except:
        print pdbname, " fail to retrive the data, need to be downloaded by hand."
        continue
    pdbfile = response.read()
    with open(completename, 'w') as f:
        f.write(pdbfile)
"""

"""
Obtain the protein sequence from pdb
"""
print "analyze pdb file"
parser = PDBParser(PERMISSIVE=1)
nameList = []
chainIDList = []
sequenceList = []
ReferenceBook = []
count = 1
# Looping for all pdb file
for filename in os.listdir("./PDB/"):
    if not filename.endswith(".pdb"):
        continue
    print "count: ", count, "pdbname: ", filename[0:4]
    count += 1
    completename = "./PDB/" + filename
    structureID = filename[0:4]
    # read pdb file
    structure = parser.get_structure(structureID, completename)
    model = structure[0]
    # Looping all the chains until find a protein chain
    for chain in list(model):
        residues = chain.get_residues()
        firstresidue = str(list(residues)[0])
        firstresname = firstresidue[9:12]
        if firstresname not in pdb2fasta:
            continue
        else:
            break
    # Extract residues and the corresponding id into list
    # residueList: Store the fasta 1-letter residue name in order
    # residueReferenceBook: Store the residue name and id together for searching in the future
    chainID = str(chain)[-2]
    chainList = chain.get_list()
    residueList = []
    residueReferenceBook = []
    for item in chainList:
        string = str(item)
        splitStr = string.split(" ")
        residueName = splitStr[1]
        if residueName not in pdb2fasta:
            continue
        residueIndex = splitStr[4]
        residueIndex = residueIndex.split("=")[1]
        residueList.append(pdb2fasta[residueName])
        residueReference = residueName + residueIndex
        residueReferenceBook.append(residueReference)

    # Collect the data into list
    nameList.append(structureID)
    chainIDList.append(chainID)
    sequenceList.append(residueList)
    ReferenceBook.append(residueReferenceBook)


"""
Search all the txt file to get the output file
"""
print "Analyze the interaction file to get the output label"
count = 1
tagSequenceList = []
for i in range(len(nameList)):
    # Prepare the necessary information for a protein sequence
    pdbname = nameList[i]
    chainID = chainIDList[i]
    residueList = sequenceList[i]
    residueReferenceBook = ReferenceBook[i]
    tagSequence = [0] * len(residueList)
    print "count: ", count, "pdbname: ", pdbname
    count += 1
    completename = "./interaction/" + pdbname + ".txt"
    with open(completename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("#"):
            continue
        string = line.split("\t")
        tmpString = string[1].split(":")
        SearchResidue = tmpString[0]
        SearchChainID = tmpString[1]
        if not SearchChainID.startswith(chainID):
            continue
        index = residueReferenceBook.index(SearchResidue)
        tagSequence[index] = 1
        res1 = residueList[index]
        res2 = pdb2fasta[SearchResidue[0:3]]
        if res1 != res2:
            print SearchResidue
    tagSequenceList.append(tagSequence)

"""
Print the result to form the data set for input and output

print "export the output and input into files "
pathname = "./dataset/"
with open(pathname+"input.txt", 'w') as f:
    for res_list in sequenceList:
        string = ' '.join(res_list)
        string = string + "\n"
        f.write(string)

with open(pathname+"output.txt", 'w') as f:
    for labels in tagSequenceList:
        labels = map(str, labels)
        string = ' '.join(labels)
        string = string + "\n"
        f.write(string)
"""