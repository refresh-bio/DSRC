#!/usr/bin/python

import pydsrc
import sys


def compress_file_illumina_lossy(infilename_, outfilename_):
    # open existing FASTQ file
    fqfile = pydsrc.FastqFile()
    fqfile.Open(infilename_)

    # create and configure DSRC archive file
    archive = pydsrc.DsrcArchive()
    archive.LossyCompression = True
    archive.TagFieldFilterMask = pydsrc.FieldMask().AddField(1).AddField(2).GetMask()
    archive.DnaCompressionLevel = 2
    archive.QualityCompressionLevel = 2
    archive.PlusRepetition = False
    archive.FastqBufferSizeMB = 256

    archive.StartCompress(outfilename_)

    # read all records from FASTQ file and write to DSRC archive
    rc = 0
    rec = pydsrc.FastqRecord()
    while fqfile.ReadNextRecord(rec):
        archive.WriteNextRecord(rec)
        rc += 1

    archive.FinishCompress()

    fqfile.Close()

    print "Success!\nRecords written: %d" % rc


def compress_file_solid_lossless(infilename_, outfilename_):
    # open existing FASTQ file
    fqfile = pydsrc.FastqFile()
    fqfile.Open(infilename_)

    # create and configure output DSRC archive file
    archive = pydsrc.DsrcArchive()
    archive.ColorSpace = True
    archive.DnaCompressionLevel = 1
    archive.QualityCompressionLevel = 1
    archive.PlusRepetition = False
    archive.FastqBufferSizeMB = 64

    archive.StartCompress(outfilename_)

    # read all records from FASTQ file and write to DSRC archive
    rc = 0
    rec = pydsrc.FastqRecord()
    while fqfile.ReadNextRecord(rec):
        archive.WriteNextRecord(rec)
        rc += 1

    archive.FinishCompress()

    fqfile.Close()

    print "Success!\nRecords written: %d" % rc


def decompress_file(infilename_, outfilename_):
    # open existing DSRC archive file
    archive = pydsrc.DsrcArchive()
    archive.StartDecompress(infilename_)

    # create results FASTQ file
    fqfile = pydsrc.FastqFile()
    fqfile.Open(outfilename_)

    # read all records from DSRC archive and write to FASTQ file
    rc = 0
    rec = pydsrc.FastqRecord()
    while archive.ReadNextRecord(rec):
        fqfile.WriteNextRecord(rec)
        rc += 1

    archive.FinishDecompress()

    fqfile.Close()

    print "Success!\nRecords written: %d" % rc


if __name__ == "__main__":
    if not (len(sys.argv) == 4 or (len(sys.argv) == 5 and sys.argv[2] == "-S")) \
       or not (sys.argv[1] == "c" or sys.argv[1] == "d"):
        print "usage: example2 <c|d> [-S] <input file> <output file>"
        exit(-1)

    try:
        if sys.argv[1] == "c":
            if len(sys.argv) == 4:
                compress_file_illumina_lossy(sys.argv[-2], sys.argv[-1])
            else:
                compress_file_solid_lossless(sys.argv[-2], sys.argv[-1])
        else:
            decompress_file(sys.argv[-2], sys.argv[-1])
    except Exception as e:
        print e.message
        exit(-1)

    exit(0)