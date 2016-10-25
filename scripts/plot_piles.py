import pandas
import matplotlib.pyplot as plt
import os, sys

for i in range(int(sys.argv[1])):
	if (os.path.isfile(str(i))):
		test = pandas.read_csv(str(i), delimiter=" ", header=None)
		test.sort_values(by=0, inplace=True)
		for j in range(len(test)):
			plt.plot((test.iloc[j,0], test.iloc[j,1]), (j,j), 'k')
		plt.savefig(str(i) + "_pile.png")
    	plt.clf()
