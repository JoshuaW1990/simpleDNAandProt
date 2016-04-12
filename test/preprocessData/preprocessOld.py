import Bio
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Selection
import csv
import urllib2
import os



"""
collect the pdb id
"""
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
Set up a dictionary for converting the pdb 3-letters residue name to fasta 1-letter residue name
"""
pdb2fasta = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}



"""
Download fasta file

print "Download fasta file"
pathname = "./fasta/"
count = 1
originalurl = "http://npidb.belozersky.msu.ru/data/pdb_new/fasta/1qna.pdb1.pdb.fas"
for pdbname in pdbidList:
    print "count: ", count, "pdbname: ", pdbname
    count += 1
    filename = pdbname + ".fas"
    completename = pathname + filename
    targeturl = originalurl.replace("1qna", pdbname)
    try:
        response = urllib2.urlopen(targeturl)
    except:
        print pdbname, " fail to retrive the data."
        continue
    fasta = response.read()
    with open(completename, 'w') as f:
        f.write(fasta)
"""

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
        print pdbname, " fail to retrive the data."
        continue
    pdbfile = response.read()
    with open(completename, 'w') as f:
        f.write(pdbfile)
"""

"""
Obtain the protein sequence from fasta

print "analyze fasta file"
nameList = []
descriptionList = []
sequenceList = []
for filename in os.listdir("./fasta/"):
    if not filename.endswith(".fas"):
        continue
    completename = "./fasta/" + filename
    readFasta = SeqIO.parse(completename, "fasta")
    for item in readFasta:
        description = item.description
        description = description.split(" ")
        if description[1] != "protein":
            continue
        else:
            chainID = description[0][-1]
            sequence = list(item.seq)
            break
    pdbname = filename[0:4]
    nameList.append(pdbname)
    descriptionList.append(chainID)
    sequenceList.append(sequence)
"""

"""
Obtain the protein sequence from pdb
"""
print "analyze pdb file"
parser = PDBParser(PERMISSIVE=1)
nameList = []
descriptionList = []
sequenceList = []
startResIndexList = []
count = 1
for filename in os.listdir("./PDB/"):
    if not filename.endswith(".pdb"):
        continue
    print "count: ", count, "pdbname: ", filename[0:4]
    count += 1
    completename = "./PDB/" + filename

    structureID = filename[0:4]
    structure = parser.get_structure(structureID, completename)
    model = structure[0]
    chainList = list(model)
    res_list = []
    for chain in chainList:
        chainID = str(chain)
        residues = chain.get_residues()
        residues = list(residues)
        # Skip this loop if this chain is a DNA chain
        firstresidue = str(residues[0])
        firstresname = firstresidue[9:12]
        if firstresname not in pdb2fasta:
            continue
        # Extract the residue sequence and print the id of chain
        startResidue = firstresidue.split("=")
        startResIndex = startResidue[2].split(" ")[0]
        startResIndex = int(startResIndex)
        for residue in residues:
            residue = str(residue)
            residueName = residue[9:12]
            if residueName not in pdb2fasta:
                break
            tmpResidue = pdb2fasta[residueName]
            res_list.append(tmpResidue)
        break
    nameList.append(structureID)
    descriptionList.append(chainID)
    startResIndexList.append(startResIndex)
    sequenceList.append(res_list)












"""
Download the corresponding interaction file

print "Download intearction file"
pathname = "./interaction/"
count = 1
originalurl = "http://npidb.belozersky.msu.ru/data/pdb_new/Hbond/1qna.pdb1.pdb.hb.txt"
for pdbname in nameList:
    print "count: ", count, "pdbname: ", pdbname
    count += 1
    filename = pdbname + ".txt"
    completename = pathname + filename
    targeturl = originalurl.replace("1qna", pdbname)
    try:
        response = urllib2.urlopen(targeturl)
    except:
        print pdbname, " fail to retrive the data."
        continue
    with open(completename, 'w') as f:
        for line in response:
            f.write(line)
"""









"""
Extract the interaction information from the txt file and return a sequence of tags for all residues in input
"""
count = 1
tagSequenceList = []
for i in range(len(nameList)):
    # Prepare the necessary information for a protein sequence
    pdbname = nameList[i]
    chainID = descriptionList[i][-2]
    sequence = sequenceList[i]
    tagSequence = [0] * len(sequence)
    startResIndex = startResIndexList[i]

    # Process to obtain the output sequence
    print "count: ", count, "pdbname: ", pdbname
    count += 1
    filename = pdbname + ".txt"
    completename = "./interaction/" + filename
    data = []
    with open(completename, 'r') as f:
        lines = f.readlines()

    # Extract the corresponding protein sequence binding information
    for line in lines:
        if line[0] == "#":
            continue
        strList = line.split("\t")
        strList = strList[1].split(":")
        if strList[1].startswith(chainID):
            string = strList[0]
            splitStr = (string[0:3], string[3:])
            data.append(splitStr)

    # Reset the tag sequence
    for item in data:
        residue = item[0]
        index = int(item[1])
        index = index - startResIndex
        tagSequence[index] = 1
        residueName = pdb2fasta[residue]
        if residueName != sequence[index]:
            print "errors in sequence checking: ", pdbname, ", incorrect index is ", (index + startResIndex)," and ", index, ", chainID is ", chainID
            print "In fasta, residue is ", sequence[index], ", in txt file, residue is ", residueName
    tagSequenceList.append(tagSequence)
