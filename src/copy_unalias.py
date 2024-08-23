#!/usr/bin/env python3
import sys
import json
import crumbs

alias_key = crumbs.get_alias_key('data/lineages.yml')

with open('/dev/stdout', 'w') as outf:
    for line in sys.stdin:
        data = line.split(',')
        lin = data[-1].strip()
        lincrumbs = crumbs.crumbs(lin, alias_key)
        outf.write(', '.join(data[:-1]) + ', ' + lin + ', True\n')
        for c in lincrumbs: outf.write(', '.join(data[:-1]) + ', ' + c + ', False\n')

