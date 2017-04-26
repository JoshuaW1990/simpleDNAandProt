"""
This method is following the article:
Identification of interface residues involved in protein-protein and protein-DNA interactions from sequence using machine learning approaches
Written by Yan, Changhui
http://lib.dr.iastate.edu/cgi/viewcontent.cgi?article=2782&context=rtd


Naive Bayes model without considering electrostatic potential
"""

from collections import defaultdict
import numpy as np
from math import sqrt


"""
Set up a dictionary for converting the pdb 3-letters residue name to fasta 1-letter residue name
"""
pdb2fasta = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


def PreprocessData():
    """
    Read the dataset from dataset directory
    input.txt: It is the multiple lines formed by residue names and each line represent a specific protein
    output.txt: It is the multiple lines formed by 1 or 0, and each line represent a specific label for the corresponding residue on the protein
    """
    pathname = "dataset/"

    # read the input.txt file
    with open(pathname + "input.txt", 'r') as f:
        lines = f.readlines()

    input_set = []
    for line in lines:
        if line.endswith("\n"):
            line = line[:-1]
        string = line.split(" ")
        startLine = ["<s>"] * 4
        endLine = ["</s>"] * 4
        string = startLine + string + endLine
        input_set.append(string)

    # read the output.txt file
    with open(pathname + "output.txt", 'r') as f:
        lines = f.readlines()

    output_set = []
    for line in lines:
        if line.endswith("\n"):
            line = line[:-1]
        string = line.split(" ")
        startLine = ['-1'] * 4
        endLine = ['-1'] * 4
        string = startLine + string + endLine
        output_set.append(string)

    # Combine the data in input file and the output file together
    combine_set = []
    for i in range(len(input_set)):
        inputList = input_set[i]
        outputList = output_set[i]
        tupleSet = []
        for j in range(len(inputList)):
            item = (inputList[j], outputList[j])
            tupleSet.append(item)
        combine_set.append(tupleSet)

    return combine_set


def BuildDataset(combine_set, index):
    """
    Build the dataset such that each line of the dataset represent an instance
    Build the inputset which is very similar to dataset except that it doesn't have the output
    """
    tmp_set = combine_set[:index] + combine_set[(index + 1):]
    test_set = combine_set[index]
    training_set = [item for sublist in tmp_set for item in sublist]
    dataset = []
    inputset = []
    outputset = []

    # loop the residue sequence in training set
    for i in range(len(training_set)):
        # each instance represent a specific sample including all the input features for feeding model
        item = training_set[i]
        if item[1] == '-1':
            continue
        instance = []
        tag = item[1]
        outputset.append(tag)
        for j in range(9):
            index = i - 4 + j
            res_name = training_set[index][0]
            instance.append(res_name)
        inputset.append(list(instance))
        instance.append(tag)
        dataset.append(list(instance))
    return (dataset, inputset, outputset, test_set)


def BuildDict(dataset):
    """
    Preprocess the dataset to obtain the counts for different features with different tags
    Store them in a dictionary
    """
    tag_dict = defaultdict(float)
    feature_dict = defaultdict(float)
    for instance in dataset:
        tag = instance[9]
        features = instance[:9]
        tag_dict[tag] += 1.0
        for i in range(len(features)):
            feature = features[i]
            item = (tag, feature, i)
            feature_dict[item] += 1.0
    return (tag_dict, feature_dict)


def Classification(theta, tag_dict, feature_dict, features):
    """
    Predict the tag with the value of theta for naive bayes model

    theta value is used to predict the tag
    """
    total = sum(tag_dict.values())
    prob = []
    for tag in range(2):
        tag = str(tag)
        probability = np.log2(tag_dict[tag] / total)
        searchFlag = 1
        for featureID in range(len(features)):
            feature = features[featureID]
            item = (tag, feature, featureID)
            if item not in feature_dict:
                searchFlag = 0
                break
            tmpProb = np.log2(feature_dict[item] / tag_dict[tag])
            probability = probability + tmpProb
        prob.append(probability)
        if searchFlag == 0:
            break
    if searchFlag == 0:
        if tag == '1':
            predict_tag = '0'
        else:
            predict_tag = '1'
    else:
        ratio = 2.0 ** (prob[1] - prob[0])
        if ratio >= theta:
            predict_tag = '1'
        else:
            predict_tag = '0'
    return predict_tag


def PredictTag(input_set, tag_dict, feature_dict, theta = 0.01):
    """
    Predict the tag with the given theta value
    """
    total = float(len(input_set))
    pred_output = []
    for features in input_set:
        # Calculate the probability for p(c=1) and p(c=0)
        pred = Classification(theta, tag_dict, feature_dict, features)
        pred_output.append(pred)
    return pred_output


def ConfusionMatrix(pred, labels):
    """
    Build the confustion matrix for assessing the preformance of the model
    """
    confusionMatrix = defaultdict(float)
    for i in range(len(labels)):
        if labels[i] == '1':
            if pred[i] == '1':
                confusionMatrix['TP'] += 1.0
            else:
                confusionMatrix['FN'] += 1.0
        else:
            if pred[i] == '1':
                confusionMatrix['FP'] += 1.0
            else:
                confusionMatrix['TN'] += 1.0
    return confusionMatrix


def CalculateCC(confusionMatrix):
    """
    calculate the correlation coefficient(CC) from the confusion matrix
    """
    TP = confusionMatrix['TP']
    TN = confusionMatrix['TN']
    FP = confusionMatrix['FP']
    FN = confusionMatrix['FN']
    cc = TP * TN - FP * FN
    tmp = sqrt((TP + FN) * (TP + FP) * (TN + FP) * (TN + FN))
    CC = float(cc) / float(tmp)
    return CC


def CalculateAccuracy(confusionMatrix):
    """
    calculate the accuracy from the confusion matrix
    """
    TP = confusionMatrix['TP']
    TN = confusionMatrix['TN']
    FP = confusionMatrix['FP']
    FN = confusionMatrix['FN']
    total = TP + TN + FP + FN
    correct = TP + TN
    accuracy = float(correct) / float(total)
    return accuracy


def Train(tag_dict, feature_dict, input_set, output_set):
    """
    Train the naive bayes model with input to obtain the best theta value as the parameter
    """
    theta = 0.01
    maxCC = 0.0
    finalTheta = 0
    confusionMatrix = None
    # use a while loop start from 0.01 for theta until 1, and the step is 0.01.
    while theta <= 1.0:
        # Obtain the predict tags
        pred_output = PredictTag(input_set, tag_dict, feature_dict, theta)
        # Buid the confusion matrix
        confusion_matrix = ConfusionMatrix(pred_output, output_set)
        # Calculate CC
        CC = CalculateCC(confusion_matrix)
        # Compare CC, keep the final theta value with the largest CC
        if CC >= maxCC:
            maxCC = CC
            finalTheta = theta
            confusionMatrix = confusion_matrix
        theta += 0.01
    return (finalTheta, confusion_matrix)


def loo_CrossValidation(combine_set):
    """
    Use the lose-one-out cross validation method for assessing the model performance
    """
    accuracy = []
    f = open("NaiveOutput1.txt", 'w')
    for i in range(len(combine_set)):
        print i
        string = "running time: " + str(i) + '\n'
        f.write(string)
        # split the dataset as input_set and output_set and train the model
        (dataset, input_set, output_set, test_set) = BuildDataset(combine_set, i)
        (tag_dict, feature_dict) = BuildDict(dataset)
        (theta, confusionMatrix) = Train(tag_dict, feature_dict, input_set, output_set)
        # Splite the test set with input and output
        print "theta = ", theta,
        trainingAccuracy = CalculateAccuracy(confusionMatrix)
        print " accuracy = ", trainingAccuracy
        string = "theta = " + str(theta) + " accuracy = " + str(trainingAccuracy) + '\n'
        f.write(string)
        test_input = []
        test_output = []
        for item in test_set:
            test_input.append(item[0])
            if item[1] == '-1':
                continue
            test_output.append(item[1])
        test_pred = []
        # predict the output according to the test input
        for i in range(len(test_input) - 8):
            res_name = test_input[i]
            startIndex = i
            endIndex = startIndex + 9
            features = test_input[startIndex:endIndex]
            pred_tag = Classification(theta, tag_dict, feature_dict, features)
            test_pred.append(pred_tag)
        confusion_matrix = ConfusionMatrix(test_pred, test_output)
        tmpAccuracy = CalculateAccuracy(confusion_matrix)
        print "test accuracy is ", tmpAccuracy
        string = "test accuracy is " + str(tmpAccuracy) + '\n'
        f.write(string)
        accuracy.append(tmpAccuracy)
    finalAccuracy = np.mean(accuracy)
    f.close()
    return finalAccuracy


combine_set = PreprocessData()
print "final accuracy after crossvalidation is: "
print loo_CrossValidation(combine_set)
