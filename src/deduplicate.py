#!/usr/bin/env python3
import sys

with open('/dev/stdout', 'w') as outf:
    bid, bdate, bline = ('', '', '')
    for line in sys.stdin:
        id = line.split('\t')
        id, date = id[0], id[2]
        if bid != id:
            outf.write(bline)
            bid, bdate, bline = (id, date, line)
        elif bdate < date:
            bid, bdate, bline = (id, date, line)
    outf.write(bline)

