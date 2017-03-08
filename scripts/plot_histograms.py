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
        test = pandas.read_csv(file_name, delimiter=" ", header=None)
        plt.plot(test[0], test[1])
        #plt.plot(test[0], test[2])
        for i in range(test[0].size):
            if (test.get_value(i, 2) != 0):
                plt.axvline(test.get_value(i, 0), color=colors[test.get_value(i, 2)])
        plt.savefig(file_name + "_hist.png")
        plt.clf()
