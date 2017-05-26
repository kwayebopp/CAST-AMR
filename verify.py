#!/usr/bin/python
import sys
i = 0
for line in open(sys.argv[1]):
    if ' ; ' in line:
        fields = line.split(' ; ')
        greaters = [1 for x in fields if len(x.split()) > 4]
        if len(fields) == len(greaters):
            i += 1
            print line
    elif ' . ' in line[:-1]:
        fields = line[:-1].split(' . ')
        if len(fields) > 0:
            i += 1
            print line

print i
