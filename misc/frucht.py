import matplotlib.pyplot, pandas, sys, seaborn

seaborn.set()

g = pandas.read_csv(sys.argv[1], header=None)
p = pandas.read_csv(sys.argv[2], header=None)

d = {}
for i, r in p.iterrows():
    d[r[0]] = [r[1], r[2]]
    matplotlib.pyplot.plot(r[1], r[2], 'bo')

for i, r in g.iterrows():
    matplotlib.pyplot.plot([d[r[0]][0], d[r[1]][0]], [d[r[0]][1], d[r[1]][1]], '-')

matplotlib.pyplot.show()
