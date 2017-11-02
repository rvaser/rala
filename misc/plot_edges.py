import pandas
import matplotlib.pyplot as plt
import os, sys

colors = [None, 'r', 'g', 'k']
for i in range(int(sys.argv[1])):
    if (os.path.isfile("e"+str(i))):
        file_name = "e"+str(i)
    else:
        file_name = ""

    if(file_name != ""):
        test = pandas.read_csv(file_name, delimiter=" ")
        plt.plot(test.iloc[:,0], test.iloc[:,1], label=str(test.columns.values[1]))
        plt.plot(test.iloc[:,0], test.iloc[:,2], label=str(test.columns.values[2]))
        plt.plot(test.iloc[:,0], test.iloc[:,3], label=str(test.columns.values[3]))
        plt.legend(loc="best")
        plt.savefig(file_name + "_edge.png")
        plt.clf()
