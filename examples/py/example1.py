#!/usr/bin/python

import pydsrc
import sys


def compress_file(infilename_, outfilename_):
    # create and configure DSRC module
    module = pydsrc.DsrcModule()
    module.LossyCompression = True
    module.TagFieldFilterMask = pydsrc.FieldMask().AddField(1).AddField(2).GetMask()
    module.DnaCompressionLevel = 2
    module.QualityCompressionLevel = 2

    module.FastqBufferSizeMB = 256
    module.ThreadsNumber = 2

    module.Compress(infilename_, outfilename_)


def decompress_file(infilename_, outfilename_):
    module = pydsrc.DsrcModule()
    module.ThreadsNumber = 2
    module.Decompress(infilename_, outfilename_)


if __name__ == "__main__":
    if len(sys.argv) != 4 or not (sys.argv[1] == "c" or sys.argv[1] == "d"):
        print "usage: example1 <c|d> <input file> <output file>"
        exit(-1)

    try:
        if sys.argv[1] == "c":
            compress_file(sys.argv[-2], sys.argv[-1])
        else:
            decompress_file(sys.argv[-2], sys.argv[-1])
    except Exception as e:
        print e.message
        exit(-1)

    print "Success!"
    exit(0)