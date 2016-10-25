import pandas
import matplotlib.pyplot as plt
import os, sys

for i in range(int(sys.argv[1])):
    if (os.path.isfile("c"+str(i))):
        file_name = "c"+str(i)
    elif (os.path.isfile("r"+str(i))):
        file_name = "r"+str(i)
    elif (os.path.isfile("l"+str(i))):
        file_name = "l"+str(i)
    else:
        file_name = ""

    if(file_name != ""):
        test = pandas.read_csv(file_name, delimiter=" ", header=None)
        plt.plot(test[0], test[1])
	plt.plot(test[0], test[2], 'r^')
        plt.savefig(file_name + "_hist.png")
        plt.clf()
