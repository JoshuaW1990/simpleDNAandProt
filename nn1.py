"""
NN without electrostatic potential
"""


from collections import defaultdict
from sklearn import svm
import numpy as np
from math import sqrt
import tensorflow as tf


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

def BuildPSSM(input_set):
    """
    Build the PSSM matrix for generating the input feature
    """
    # Initialize the pssm
    pssm = [{} for i in range(9)]
    total_res_dict = {}
    residue_list = pdb2fasta.values() + ['<s>', '</s>']
    for residue in residue_list:
        for i in range(len(pssm)):
            total_res_dict[residue] = 0.0
            pssm[i][residue] = 0.0

    # Count the residue
    for instance in input_set:
        for i in range(len(instance)):
            residue = instance[i]
            pssm[i][residue] += 1.0
            total_res_dict[residue] += 1.0

    # Calculate the log value
    total_res_num = sum(total_res_dict.values())
    for residue in residue_list:
        total_res_dict[residue] = total_res_dict[residue] / float(total_res_num)
    for res_dict in pssm:
        for residue in residue_list:
            frequency = float(res_dict[residue] + total_res_dict[residue]) / float(len(input_set) + 1)
            res_dict[residue] = np.log(frequency / total_res_dict[residue])

    return pssm


def preprocess_data(pssm, input_set, output_set, test_set):
    """
    Preprocess the dataset such that use the value in pssm to replace the residue and form input features
    """
    X_train = []
    Y_train = []
    X_test = []
    Y_test = []
    # generate the training set
    for i in range(len(input_set)):
        instance = input_set[i]
        if output_set[i] == '0':
            tag = [1., 0.]
        else:
            tag = [0., 1.]
        Y_train.append(tag)
        input_feature = []
        for i in range(len(instance)):
            residue = instance[i]
            input_feature.append(pssm[i][residue])
        X_train.append(input_feature)

    # generate the test set
    for i in range(len(test_set)):
        item = test_set[i]
        if item[1] == '-1':
            continue
        input_feature = []
        if output_set[i] == '0':
            tag = [1., 0.]
        else:
            tag = [0., 1.]
        Y_test.append(tag)
        for j in range(9):
            index = i - 4 + j
            res_name = test_set[index][0]
            input_feature.append(pssm[j][res_name])
        X_test.append(input_feature)

    return (X_train, Y_train, X_test, Y_test)


def Cross_validation():
    """
    leave one out cross validation for the dataset
    """
    test_accuracy = []
    combine_set = PreprocessData()
    f = open("ann1.txt", 'w')
    for i in range(len(combine_set)):
        print i
        string = "running time: " + str(i) + '\n'
        f.write(string)
        # prepare the training set and the test set
        (dataset, input_set, output_set, test_set) = BuildDataset(combine_set, i)
        pssm = BuildPSSM(input_set)
        (X_train, Y_train, X_test, Y_test) = preprocess_data(pssm, input_set, output_set, test_set)

        # training the model with neural network in tensorflow
        x = tf.placeholder("float", [None, 9])
        W = tf.Variable(tf.zeros([9, 2]))
        b = tf.Variable(tf.zeros([2]))
        y = tf.nn.softmax(tf.matmul(x, W) + b)
        y_ = tf.placeholder("float", [None, 2])
        cross_entropy = -tf.reduce_sum(y_*tf.log(y))
        train_step = tf.train.GradientDescentOptimizer(0.01).minimize(cross_entropy)
        init = tf.initialize_all_variables()
        sess = tf.Session()
        sess.run(init)
        for i in range(1000):
            sess.run(train_step, feed_dict={x: X_train, y_: Y_train})

        # predict and calculate the accuracy
        correct_prediction = tf.equal(tf.argmax(y,1), tf.argmax(y_,1))
        accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))
        test_acc = sess.run(accuracy, feed_dict={x: X_test, y_: Y_test})
        sess.close()
        print "test accuracy: ", test_acc
        string = "test accuracy is " + str(test_acc) + '\n'
        f.write(string)
        test_accuracy.append(test_acc)
    print np.mean(test_accuracy)
    f.close()
    return test_accuracy


test_accuracy = Cross_validation()





