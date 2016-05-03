import numpy as np
import scipy.stats as stat
import csv

with open("NaiiveOutput.txt", 'r') as f:
    lines = f.readlines()


time = []
test_acc = []
train_acc = []
theta = []
for line in lines:
    if line.startswith("running"):
        string = line.split(" ")
        time.append(int(string[-1].split("\n")[0]))
        continue
    if line.startswith("theta"):
        string = line.split("=")
        theta.append(float(string[1].split(" ")[1]))
        train_acc.append(float(string[-1].split("\n")[0]))
        continue
    if line.startswith("test"):
        string = line.split(" ")
        test_acc.append(float(string[-1].split("\n")[0]))
        continue
    else:
        continue

time = np.array(time)
test_acc = np.array(test_acc)
train_acc = np.array(train_acc)
theta = np.array(theta)

with open("data.csv", "w") as csvfile:
    write = csv.writer(csvfile, dialect = 'excel')
    for i in range(len(time)):
        write.writerow([time[i], test_acc[i], train_acc[i], theta[i]])