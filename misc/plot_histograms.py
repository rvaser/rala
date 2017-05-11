import pandas
import matplotlib.pyplot as plt
import os, sys

colors = [None, 'r', 'g', 'k']
for i in range(int(sys.argv[1])):
    if (os.path.isfile("c"+str(i))):
        file_name = "c"+str(i)
    elif (os.path.isfile("h"+str(i))):
        file_name = "h"+str(i)
    else:
        file_name = ""

    if(file_name != ""):
        test = pandas.read_csv(file_name, delimiter=" ")
        plt.plot(test.iloc[:,0], test.iloc[:,1], label=str(test.columns.values[1]))
        plt.plot(test.iloc[:,0], test.iloc[:,3], label=str(test.columns.values[3]))
        plt.plot(test.iloc[:,0], test.iloc[:,4], label=str(test.columns.values[4]))
        for i in range(test.iloc[:,0].size):
            if (test.iloc[i, 2] != 0):
                plt.axvline(test.iloc[i,0], color=colors[test.iloc[i, 2]])
        plt.legend(loc="best")
        plt.savefig(file_name + "_hist.png")
        plt.clf()
