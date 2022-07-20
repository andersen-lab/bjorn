import sys
import zstandard
import argparse
from base64io import Base64IO

args = argparse.ArgumentParser(description='Zstandard compress a file line-wise')
args.add_argument('--level', type=int, help='compression level')
args.add_argument('--dict', type=argparse.FileType('rb'), help='dictionary to use')
args = vars(args.parse_args())

cctx = zstandard.ZstdCompressor(level=args['level'], dict_data=zstandard.ZstdCompressionDict(args['dict'].read()))
with open('/dev/stdout', 'wb') as outf, Base64IO(outf) as base64:
    comp = cctx.stream_writer(base64)
    for line in sys.stdin:
        line = line.split('\t')
        if len(line) > 1: outf.write(bytearray('\t'.join(line[:-1])+'\t', 'utf8'))
        comp.write(bytearray(line[-1], 'utf8'))
        comp.flush(zstandard.FLUSH_FRAME)
        base64.close()
        base64.closed = False
        outf.write(bytearray('\n', 'utf8'))
        outf.flush()
