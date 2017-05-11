# changes read names to integers representing their id (starting from one) and outputs the result to stdout
import sys

id = 1
i = 0
with (open(sys.argv[1])) as f:
    for line in f:
        if (i % 4 == 0):
            print("@%d" % id)
            id += 1
        else:
            print(line.rstrip())
        i += 1
