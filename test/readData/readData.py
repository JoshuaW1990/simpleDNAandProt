pathname = "../preprocessData/dataset/"

with open(pathname + "input.txt", 'r') as f:
    lines = f.readlines()

input_set = []
for line in lines:
    if line.endswith("\n"):
        line = line[:-1]
    string = line.split(" ")
    input_set.append(string)

with open(pathname + "output.txt", 'r') as f:
    lines = f.readlines()

output_set = []
for line in lines:
    if line.endswith("\n"):
        line = line[:-1]
    string = line.split(" ")
    output_set.append(string)

"""
Check the length of the inputset and output_set
"""
for i in range(len(input_set)):
    if len(input_set[i]) != len(output_set[i]):
        print "incorrect length of input and output list"
