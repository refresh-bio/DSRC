#!/usr/bin/python

import pydsrc
import sys


def compress_file(infilename_, outfilename_):
    # open existing FASTQ file
    fqfile = pydsrc.FastqFileRecordsReader()
    fqfile.Open(infilename_)

    # define compression configuration settings for DSRC module:
    # lossless mode (by default) with 'fast' performance setting
    # (corresponds to '-m0' switch in DSRC binary)
    settings = pydsrc.CompressionSettings()
    settings.DnaCompressionLevel = 0
    settings.QualityCompressionLevel = 0
    settings.FastqBufferSizeMb = 8

    # create DSRC archive
    archive = pydsrc.DsrcArchiveRecordsWriter()
    archive.StartCompress(outfilename_, settings)

    # read records from FASTQ file writing to DSRC archive
    rc = 0
    rec = pydsrc.FastqRecord()
    while fqfile.ReadNextRecord(rec):
        archive.WriteNextRecord(rec)
        rc += 1

    archive.FinishCompress()

    fqfile.Close()

    print "Success!\nRecords written: %d" % rc


def compress_file_lossy(infilename_, outfilename_):
    # open existing FASTQ file
    fqfile = pydsrc.FastqFileRecordsReader()
    fqfile.Open(infilename_)

    # define compression configuration settings for DSRC module:
    # use 'max' ratio compression mode
    # (corresponds to '-m2' switch in DSRC binary)
    settings = pydsrc.CompressionSettings()
    settings.DnaCompressionLevel = 3
    settings.QualityCompressionLevel = 2
    settings.FastqBufferSizeMb = 256

    # set lossy mode with Illumina 8-binning qualities scheme
    settings.LossyQualityCompression = True

    # keep only 2 first tokens in read tag field
    settings.TagFieldFilterMask = pydsrc.FieldMask().AddField(1).AddField(2).GetMask()

    # create DSRC archive
    archive = pydsrc.DsrcArchiveRecordsWriter()
    archive.StartCompress(outfilename_, settings)

    # read records from FASTQ file writing to DSRC archive
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
    archive = pydsrc.DsrcArchiveRecordsReader()
    archive.StartDecompress(infilename_)

    # create results FASTQ file
    fqfile = pydsrc.FastqFileRecordsWriter()
    fqfile.Open(outfilename_)

    # read all records from DSRC archive writing to FASTQ file
    rc = 0
    rec = pydsrc.FastqRecord()
    while archive.ReadNextRecord(rec):
        fqfile.WriteNextRecord(rec)
        rc += 1

    archive.FinishDecompress()

    fqfile.Close()

    print "Success!\nRecords written: %d" % rc


if __name__ == "__main__":
    if not (len(sys.argv) == 4 or (len(sys.argv) == 5 and sys.argv[2] == "L")) \
       or not (sys.argv[1] == "c" or sys.argv[1] == "d"):
        print "usage: example2 <c|d> [L] <input file> <output file>"
        exit(-1)

    lossyMode = (len(sys.argv) == 5) and sys.argv[2] == "L"

    try:
        if sys.argv[1] == "c":
            if lossyMode:
                compress_file_lossy(sys.argv[-2], sys.argv[-1])
            else:
                compress_file(sys.argv[-2], sys.argv[-1])
        else:
            decompress_file(sys.argv[-2], sys.argv[-1])

    except Exception as e:
        print e.message
        exit(-1)
