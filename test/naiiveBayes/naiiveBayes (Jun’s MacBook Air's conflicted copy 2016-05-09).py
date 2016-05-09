"""
This method is following the article:
Identification of interface residues involved in protein-protein and protein-DNA interactions from sequence using machine learning approaches
Written by Yan, Changhui
http://lib.dr.iastate.edu/cgi/viewcontent.cgi?article=2782&context=rtd
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


"""
Read the dataset
"""
def PreprocessData():
    pathname = "../preprocessData/dataset/"

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

"""
Build the dataset such that each line of the dataset represent an instance
Build the inputset which is very similar to dataset except that it doesn't have the output
"""
def BuildDataset(combine_set, index):
    tmp_set = combine_set[:index] + combine_set[(index + 1):]
    test_set = combine_set[index]
    training_set = [item for sublist in tmp_set for item in sublist]
    dataset = []
    inputset = []
    outputset = []
    for i in range(len(training_set)):
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


"""
Preprocess the dataset to obtain the counts for different features with different tags
Store them in a dictionary
"""
def BuildDict(dataset):
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

"""
Train the value of theta
Return theta
"""
def Classification(theta, tag_dict, feature_dict, features):
    predict_tag = None
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
    total = float(len(input_set))
    pred_output = []
    for features in input_set:
        # Calculate the probability for p(c=1) and p(c=0)
        pred = Classification(theta, tag_dict, feature_dict, features)
        pred_output.append(pred)
    return pred_output

def ConfusionMatrix(pred, labels):
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
    TP = confusionMatrix['TP']
    TN = confusionMatrix['TN']
    FP = confusionMatrix['FP']
    FN = confusionMatrix['FN']
    cc = TP * TN - FP * FN
    tmp = sqrt((TP + FN) * (TP + FP) * (TN + FP) * (TN + FN))
    CC = float(cc) / float(tmp)
    return CC

def CalculateAccuracy(confusionMatrix):
    TP = confusionMatrix['TP']
    TN = confusionMatrix['TN']
    FP = confusionMatrix['FP']
    FN = confusionMatrix['FN']
    total = TP + TN + FP + FN
    correct = TP + TN
    accuracy = float(correct) / float(total)
    return accuracy

def Train(tag_dict, feature_dict, input_set, output_set):
    theta = 0.01
    maxCC = 0.0
    finalTheta = 0
    confusionMatrix = None
    while theta <= 1.0:
        # Obtain the predict tags
        pred_output = PredictTag(input_set, tag_dict, feature_dict, theta)
        # Buid the confusion matrix
        confusion_matrix = ConfusionMatrix(pred_output, output_set)
        # Calculate CC
        CC = CalculateCC(confusion_matrix)
        # Compare CC
        if CC >= maxCC:
            maxCC = CC
            finalTheta = theta
            confusionMatrix = confusion_matrix
        theta += 0.01
    return (finalTheta, confusion_matrix)





def loo_CrossValidation(combine_set):
    accuracy = []
    f = open("NaiiveOutput.txt", 'w')
    for i in range(len(combine_set)):
        print i
        string = "running time: " + str(i) + '\n'
        f.write(string)
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
    f.write(str(finalAccuracy))
    f.close()
    return finalAccuracy



combine_set = PreprocessData()

(dataset, input_set, output_set, test_set) = BuildDataset(combine_set, 55)
"""
(tag_dict, feature_dict) = BuildDict(dataset)
(theta, confusionMatrix) = Train(tag_dict, feature_dict, input_set, output_set)

print "final accuracy after crossvalidation is: "
print loo_CrossValidation(combine_set)
"""
