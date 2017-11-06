import sys

glen = int(sys.argv[2])
id = 0
arr = []
with (open(sys.argv[1])) as f:
    for line in f:
        if (line[0] == '>'):
            if (id > 0):
                arr.append(len(data))
            name = line.rstrip()
            data = ""
            id += 1
        else:
            data += line.rstrip()

arr.append(len(data))

arr.sort(reverse=True)
sum = 0
for i in range(len(arr)):
    sum += arr[i]
    if sum > glen / 2:
        print(arr[i])
        break
