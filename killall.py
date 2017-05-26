#!/usr/bin/python
import sys
import os
nodes = sys.argv[1].split(',')
for node in nodes:
    os.system(r'ssh %s "killall -u xpeng"' % node)
