#!/usr/bin/env python3
import sys
import json

with open('../data/alias_key.json', 'r') as alias_key:
    alias_key = json.load(alias_key)
alias_key = [(n, {v:k for k,v in alias_key.items() if isinstance(v, str) and len(v.split('.'))==n}) for n in range(16)]
alias_key = list(reversed([(n,d) for n,d in alias_key if len(d) > 0]))

with open('/dev/stdout', 'w') as outf:
    for line in sys.stdin:
        head = line.split('.')
        tail = [head.pop()]
        for n,key in alias_key:
            if n > len(head): continue
            while len(head) > n:
                tail = [head.pop()] + tail
            if '.'.join(head) in key:
                head = [key['.'.join(head)]]
                break
        outf.write('.'.join(head+tail))
