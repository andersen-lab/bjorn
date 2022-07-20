import sys
import zstandard
import argparse
import base64

args = argparse.ArgumentParser(description='Zstandard compress a file line-wise')
args.add_argument('--dict', type=argparse.FileType('rb'), help='dictionary to use')
args = vars(args.parse_args())

cctx = zstandard.ZstdDecompressor(dict_data=zstandard.ZstdCompressionDict(args['dict'].read()))
with open('/dev/stdout', 'wb') as outf, cctx.stream_writer(outf) as z:
    for line in sys.stdin:
        line = line.split('\t')
        if len(line) > 1: outf.write(bytearray('\t'.join(line[:-1])+'\t', 'utf8'))
        z.write(base64.b64decode(line[-1]))
        z.flush()
        outf.flush()
