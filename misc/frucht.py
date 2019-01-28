import matplotlib.pyplot, pandas, sys, seaborn

seaborn.set()
seaborn.set_style("white")
scpg = seaborn.cubehelix_palette(rot=-.4)
scpr = seaborn.color_palette("Reds")

p = pandas.read_csv(sys.argv[1], header=None)
e = pandas.read_csv(sys.argv[2], header=None)
t = pandas.read_csv(sys.argv[3], header=None)

matplotlib.pyplot.figure(figsize=(16,16))

d = {}
for i, r in p.iterrows():
    d[r[0]] = [r[1], r[2], r[3]]

for i, r in e.iterrows():
    clr = scpg[2]
    if (d[r[0]][2] == 1 or d[r[1]][2] == 1):
        clr = scpr[3]
    matplotlib.pyplot.plot([d[r[0]][0], d[r[1]][0]], [d[r[0]][1], d[r[1]][1]], '-', c=clr)

for i, r in t.iterrows():
    matplotlib.pyplot.plot([d[r[0]][0], d[r[1]][0]], [d[r[0]][1], d[r[1]][1]], ':', c=scpg[1])

for i, r in p.iterrows():
    matplotlib.pyplot.plot(r[1], r[2], '.', c=scpg[4], markersize=(5 if r[4] == 1 else 25))

matplotlib.pyplot.xticks([])
matplotlib.pyplot.yticks([])
matplotlib.pyplot.savefig('graph.svg', format='svg', dpi=1200)
