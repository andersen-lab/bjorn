#!/usr/bin/env python3
import sys
import json

with open('alias_key.json', 'r') as alias_key:
    alias_key = json.load(alias_key)
alias_key = {k:v for k,v in alias_key.items() if isinstance(v, str)}
#alias_key = [(n, {v:k for k,v in alias_key.items() if isinstance(v, str) and len(v.split('.'))==n}) for n in range(16)]
#alias_key = list(reversed([(n,d) for n,d in alias_key if len(d) > 0]))

with open('/dev/stdout', 'w') as outf:
    for line in sys.stdin:
        leaf = True
        data = line.split(',')
        head = data[-1].strip()
        while len(head) > 0:
            outf.write(', '.join(data[:-1]) + ', ' + head + ', ' + str(leaf) + '\n')
            if head in alias_key:
                head = alias_key[head]
#            else:
            head = '.'.join(head.split('.')[:-1])
            leaf = False
