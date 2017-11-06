# changes read names to integers representing their id (starting from one) and outputs the result to stdout
# usage: python fasta_formatter.py <input FASTA file>

import sys

id = 0
with (open(sys.argv[1])) as f:
    for line in f:
        if (line[0] == '>'):
            if (id > 0):
                print(">%d" % id)
                print(data)
            data = ""
            id += 1
        else:
            data += line.rstrip()

print(">%d" % id)
print(data)
