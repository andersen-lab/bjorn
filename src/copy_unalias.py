#!/usr/bin/env python3
import sys
import json
import crumbs

#with open('alias_key.json', 'r') as alias_key:
#    alias_key = json.load(alias_key)
#alias_key = {k:v for k,v in alias_key.items() if isinstance(v, str)}
#alias_key = [(n, {v:k for k,v in alias_key.items() if isinstance(v, str) and len(v.split('.'))==n}) for n in range(16)]
#alias_key = list(reversed([(n,d) for n,d in alias_key if len(d) > 0]))

alias_key = crumbs.get_alias_key('/data/lineages.yml')

with open('/dev/stdout', 'w') as outf:
    for line in sys.stdin:
        data = line.split(',')
        lin = data[-1].strip()
        lincrumbs = crumbs.crumbs(lin, alias_key)
        outf.write(', '.join(data[:-1]) + ', ' + lin + ', True\n')
        for c in lincrumbs: outf.write(', '.join(data[:-1]) + ', ' + c + ', False\n')

