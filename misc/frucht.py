import matplotlib.pyplot, pandas, sys, seaborn

seaborn.set()
seaborn.set_style("white")
scp = seaborn.cubehelix_palette(rot=-.4)

p = pandas.read_csv(sys.argv[1], header=None)
e = pandas.read_csv(sys.argv[2], header=None)
t = pandas.read_csv(sys.argv[3], header=None)

d = {}
for i, r in p.iterrows():
    d[r[0]] = [r[1], r[2]]

for i, r in e.iterrows():
    matplotlib.pyplot.plot([d[r[0]][0], d[r[1]][0]], [d[r[0]][1], d[r[1]][1]], '-', c=scp[2])

for i, r in t.iterrows():
    matplotlib.pyplot.plot([d[r[0]][0], d[r[1]][0]], [d[r[0]][1], d[r[1]][1]], ':', c=scp[1])

for i, r in p.iterrows():
    matplotlib.pyplot.plot(r[1], r[2], '.', c=scp[4])

matplotlib.pyplot.xticks([])
matplotlib.pyplot.yticks([])
matplotlib.pyplot.show()
