#!/usr/bin/python

import os
import sys

# set this up to a local DSRC binary
dsrc_exec = "./dsrc"


# command patterns to be executed
dsrc_compress_cmd = "%s c {params} {infile} {outfile}" % dsrc_exec
dsrc_decompress_cmd = "%s d {params} {infile} {outfile}" % dsrc_exec

dsrc_compress_stdio_cmd = "%s c -s {params} {outfile} < {infile}" % dsrc_exec
dsrc_decompress_stdio_cmd = "%s d -s {params} {infile} > {outfile}" % dsrc_exec

diff_cmd = "diff -q {infile} {outfile}"


class RunException(Exception):
    def __init__(self, msg_):
         self.msg = msg_
    def __str__(self):
        return repr(self.msg)

def run_cmd(cmd_):
    # TODO: this one should be done using subprocess module, 
    # but handling properly both stdout/stdin and stderr is a mess...
    ret = os.system(cmd_)
    if ret != 0:
        raise RunException("Error running command: %s (exit status: $d)" % (cmd_, ret))

def perform_test(params_, infile_, checkOutput_ = True, usesStdio_ = False):
    dsrc_tmp_file = "__out.dsrc"
    fastq_tmp_file = "__out.fastq"
    
    # set-up : cleanup
    if os.path.isfile(dsrc_tmp_file):
        os.remove(dsrc_tmp_file)
        
    if os.path.isfile(fastq_tmp_file):
        os.remove(fastq_tmp_file)
    
    # perform test
    try:
        print "Compressing..."
        if usesStdio_:
            pattern = dsrc_compress_stdio_cmd
        else:
            pattern = dsrc_compress_cmd
        
        cmd = pattern.format(params=params_,
                             infile=infile_,
                             outfile=dsrc_tmp_file)
        run_cmd(cmd)
        
        print "Decompressing..."
        if usesStdio_:
            pattern = dsrc_decompress_stdio_cmd
        else:
            pattern = dsrc_decompress_cmd
            
        cmd = pattern.format(params=params_,
                             infile=dsrc_tmp_file,
                             outfile=fastq_tmp_file)
        run_cmd(cmd)
        
        if checkOutput_:
            print "File check..."
            cmd = diff_cmd.format(infile=infile_,
                                  outfile=fastq_tmp_file)
            run_cmd(cmd)
    
    except RunException as exc:
        print "FAIL: " + str(exc)
    else:
        print "PASS"
    
    # tear-down : cleanup
    #
    if os.path.isfile(dsrc_tmp_file):
        os.remove(dsrc_tmp_file)

    if os.path.isfile(fastq_tmp_file):
        os.remove(fastq_tmp_file)
    

def run_tests(infile_):
    
    if not os.path.isfile(infile_):
        print "Error: input file does not exist: " + infile_
        return
    
    print "Running tests on %s file..." % infile_
    
    # simple test : only lossless mode
    #
    for th in ["-t1", "-t4"]:
        for m in ["-m0", "-m1", "-m2"]:
            params = "%s %s" % (th, m)
            print "**** Running case: (%s + file+file) ****" % params 
            perform_test(params, infile_, True)
            print "**** Running case: (%s + file+stdio) ****" % params 
            perform_test(params, infile_, True, True)
            
    # TODO: exhaustive test
    # * different compression modes + buffer sizes
    # * lossy quality mode
    # * lossy read_id mode
    # * performance tests


if __name__ == "__main__":
    if len(sys.argv) != 2 or "-h" in sys.argv[1]:
        print "usage: test_dsrc_bin.py <FASTQ filename>"
        exit(0)
        
    run_tests(sys.argv[1])
