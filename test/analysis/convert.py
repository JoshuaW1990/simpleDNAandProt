import numpy as np
import scipy.stats as stats
import csv

accuracy = []
# Read Naive Bayes model
for i in range(1, 3):
    filename = "NaiveOutput" + str(i) + ".txt"
    test_accuracy = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("test"):
            string = line.split(" ")
            test_acc = string[-1][:-1]
            test_acc = float(test_acc)
            test_accuracy.append(test_acc)
    accuracy.append(test_accuracy)

# Read svm
for i in range(1, 3):
    filename = "svm" + str(i) + ".txt"
    test_accuracy = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("test"):
            string = line.split(" ")
            test_acc = string[-1][:-1]
            test_acc = float(test_acc)
            test_accuracy.append(test_acc)
    accuracy.append(test_accuracy)

# Read Neural Network
for i in range(1, 3):
    filename = "ann" + str(i) + ".txt"
    test_accuracy = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("test"):
            string = line.split(" ")
            test_acc = string[-1][:-1]
            test_acc = float(test_acc)
            test_accuracy.append(test_acc)
    accuracy.append(test_accuracy)


(tValue, pValue) = stats.ttest_ind(accuracy[1], accuracy[3])
print tValue
print pValue
mean_accuracy = np.mean(accuracy, axis=1)
print "finish the program"