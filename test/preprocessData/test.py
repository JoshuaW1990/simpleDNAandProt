from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Selection

"""
Set up a dictionary for converting the pdb 3-letters residue name to fasta 1-letter residue name
"""
pdb2fasta = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


parser = PDBParser(PERMISSIVE=1)
structureID = "1a36"
completename = "./PDB/" + "1a36.pdb"
structure = parser.get_structure(structureID, completename)
model = structure[0]
for chain in list(model):
    residues = chain.get_residues()
    firstresidue = str(list(residues)[0])
    firstresname = firstresidue[9:12]
    if firstresname not in pdb2fasta:
        continue
    else:
        break

chainID = str(chain)[-2]
chainList = chain.get_list()
residueList = []
residueReferenceBook = []
for item in chainList:
    string = str(item)
    splitStr = string.split(" ")
    residueName = splitStr[1]
    residueIndex = splitStr[4]
    residueIndex = residueIndex.split("=")[1]
    if residueName not in pdb2fasta:
        continue
    residueList.append(pdb2fasta[residueName])
    residueReference = residueName + residueIndex
    residueReferenceBook.append(residueReference)

"""
Test txt file
"""
tagSequenceList = [0] * len(residueList)
completename = "./interaction/" + "1a36.txt"
with open(completename, 'r') as f:
    lines = f.readlines()
for line in lines:
    if line.startswith("#"):
        continue
    string = line.split("\t")
    SearchResidue = string[1].split(":")[0]
    index = residueReferenceBook.index(SearchResidue)
    tagSequenceList[index] = 1
    res1 = residueList[index]
    res2 = pdb2fasta[SearchResidue[0:3]]
    if res1 != res2:
        print SearchResidue
