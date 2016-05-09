"""
This method is following the article:
Identification of interface residues involved in protein-protein and protein-DNA interactions from sequence using machine learning approaches
Written by Yan, Changhui
http://lib.dr.iastate.edu/cgi/viewcontent.cgi?article=2782&context=rtd
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

# Get the electrostatic potential according to the residue name
values = pdb2fasta.values()
fasta2elec = dict.fromkeys(values, 0.5)
fasta2elec[pdb2fasta['ARG']] = 1.0
fasta2elec[pdb2fasta['LYS']] = 1.0
fasta2elec[pdb2fasta['HIS']] = 1.0
fasta2elec[pdb2fasta['ASP']] = 0.0
fasta2elec[pdb2fasta['GLU']] = 0.0
fasta2elec['<s>'] = 0.5
fasta2elec['</s>'] = 0.5

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

def BuildPSSM(input_set):
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

def generate_Instance_Vector(pssm, instance):
    instance_vector = []
    for i in range(len(instance)):
        residue = instance[i]
        instance_vector.append(pssm[i][residue])
    return instance_vector

def preprocess_data(pssm, input_set, output_set, test_set):
    X_train = []
    Y_train = []
    X_test = []
    Y_test = []
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
            input_feature.append(fasta2elec[residue])
        X_train.append(input_feature)

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
            input_feature.append(fasta2elec[res_name])
        X_test.append(input_feature)

    return (X_train, Y_train, X_test, Y_test)

"""helper function
"""
#accuracy
def accurayc(pred_labels, labels):
    correct = 0
    total = len(labels)
    for i in range(len(labels)):
        pred = pred_labels[i]
        label = labels[i]
        if pred == label:
            correct += 1
    return float(correct) / float(total)

"""leave one out cross validation for the dataset
"""
def Cross_validation():
    test_accuracy = []
    combine_set = PreprocessData()
    f = open("ann2.txt", 'w')
    for i in range(len(combine_set)):
        print i
        string = "running time: " + str(i) + '\n'
        f.write(string)
        (dataset, input_set, output_set, test_set) = BuildDataset(combine_set, i)
        pssm = BuildPSSM(input_set)
        (X_train, Y_train, X_test, Y_test) = preprocess_data(pssm, input_set, output_set, test_set)
        x = tf.placeholder("float", [None, len(X_train[0])])
        W = tf.Variable(tf.zeros([len(X_train[0]), 2]))
        b = tf.Variable(tf.zeros([2]))
        y = tf.nn.softmax(tf.matmul(x,W) + b)

        y_ = tf.placeholder("float", [None,2])
        cross_entropy = -tf.reduce_sum(y_*tf.log(y))

        train_step = tf.train.GradientDescentOptimizer(0.01).minimize(cross_entropy)

        init = tf.initialize_all_variables()
        sess = tf.Session()
        sess.run(init)

        for i in range(1000):
            sess.run(train_step, feed_dict={x: X_train, y_: Y_train})

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





"""Using tensorflow
"""

"""Base model

x = tf.placeholder("float", [None, 9])
W = tf.Variable(tf.zeros([9, 2]))
b = tf.Variable(tf.zeros([2]))
y = tf.nn.softmax(tf.matmul(x,W) + b)

y_ = tf.placeholder("float", [None,2])
cross_entropy = -tf.reduce_sum(y_*tf.log(y))

train_step = tf.train.GradientDescentOptimizer(0.01).minimize(cross_entropy)

init = tf.initialize_all_variables()
sess = tf.Session()
sess.run(init)

for i in range(1000):
    if i % 100 == 0:
        print i
    sess.run(train_step, feed_dict={x: X_train, y_: Y_train})



correct_prediction = tf.equal(tf.argmax(y,1), tf.argmax(y_,1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))
acc = sess.run(accuracy, feed_dict={x: X_test, y_: Y_test})

sess.close()

print acc


test accuracy:
0.886029
"""




