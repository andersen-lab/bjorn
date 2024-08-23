#!/usr/bin/env python3
import sys
import json
import hashlib

with open('/dev/stdout', 'w') as outf:
    for line in sys.stdin:
        data = json.loads(line)
        data['id'] = hashlib.md5(data['id'].encode('utf8')).hexdigest()
        print(json.dumps(data))
