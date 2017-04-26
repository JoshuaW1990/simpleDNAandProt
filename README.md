# README


## Introduction

In this project, naive bayes model, support vector machine and neural network are implemented for recognizing the DNA binding sites on protein sequence. Through leave one out cross validation and t-test, it is indicated that support vector machine has the best performance, while the performance of the naive bayes model is the worst

## Data:

### Data Collection

The list of the proteins were obtained from the [article](http://nar.oxfordjournals.org/content/31/24/7189),  which includes 56 different types DNA-protein complexes. They are listed as below:

![image](https://github.com/JoshuaW1990/simpleDNAandProt/blob/master/pictures/proteins.PNG?raw=true)

All the pdb file and the interaction file were downloaded from [NPIDB](http://npidb.belozersky.msu.ru/).

### Data Preprocess

#### Generate the residue sequence of the protein from the pdb file

1. Extract the residue sequence from the pdb file into a list to represent a chain of proteins
2. Convert the residue name to fasta 1-letter residue name

#### Generate the label sequence of the protein from the pdb file and the interactive file

In the interactive file, it records all the residues whose distances to the DNA molecules are smaller than 3.5 Angstrom.

1. Check the interation file and the corresponding chain sequence of the specific protein
2. Filter the interaction file to remove the unrelated record
3. Mark the label of the specific residue indicated in the interaction file as 1. The other labels are 0
4. For each sequence of the proteins, add four fake residues at the beginning and add four fake residues at the end of the beginning. In this way, each real residue in the sequence has the same number of neighboring residues for feature extraction (see the figure below).

![image](https://github.com/JoshuaW1990/simpleDNAandProt/blob/master/pictures/residue2.PNG?raw=true)

#### Features

After obtaining the dataset, it is also critical to determine the features for the instance in the models. Here, two different types of the features are considered: 

1. the neighboring residue sequence. 
2. the electrostatic potential of the residues in the neighboring residues.

**Neighboring residue sequence**

The length of the neighboring residue sequence will be 9. In this neighboring residue sequence, the first four residues are the residues before the target residue and the last four residues are the residues after the target residue. The middle residue is the target residue itself(See the figure below).

![image](https://github.com/JoshuaW1990/simpleDNAandProt/blob/master/pictures/residue.PNG?raw=true)

**Electrostatic Potential**

There are many ways to calculate the electrostatic potential. A simple way is to assign the electrostatic potential to -1, 0, or 1 according to the specific type of the residues (see the table below).

| Discrete Value of electrostatic potential | residues |
| ----------------------------------------- | :------- |
| Positive | Arg, Lys, His |
| Negative | Asp, Glu |
| Neutral | All others |

We also still consider the electrostatic potential of 9 residues for each residue in the model.

**Postition-specific scoring matrix(PSSM)**

In Support Vector Machine(SVM) and Neural Network (NN), it is necessary to convert the features of each residue which are residue sequence into the float numbers. Some researchers indicated that using Postition-specific scoring matrix(PSSM) to represent the neighboring residue sequence is a promising approach. PSSM represents the probabilities of the residues appeared at a specific position(see the figure below). 

![image](https://github.com/JoshuaW1990/simpleDNAandProt/blob/master/pictures/pssm.PNG?raw=true)

After building the PSSM for the training set, we convert the neighboring residue sequence of each residue to a number sequence by this matrix. This number sequence represents a vector of features for a instance.

## Techniques

- Naive Bayes model ([Reference](http://lib.dr.iastate.edu/cgi/viewcontent.cgi?article=2782&context=rtd))
- Support Vector machine ([Reference](http://nar.oxfordjournals.org/content/37/suppl_2/W396))
- Artificial Neural Network ([Reference](http://nar.oxfordjournals.org/content/35/5/1465))




