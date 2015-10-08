import os
import sys

# set this up to a local DSRC binary
dsrc_exec = "./dsrc"
lib1_exec = "./example1"
lib2_exec = "./example2"
lib3_exec = "./example3"


# command patterns to be executed
dsrc_compress_cmd = "%s c {params} {infile} {outfile}" % dsrc_exec
dsrc_decompress_cmd = "%s d {infile} {outfile}" % dsrc_exec

lib1_compress_cmd = "%s c {params} {infile} {outfile}" % lib1_exec
lib1_decompress_cmd = "%s d {infile} {outfile}" % lib1_exec

lib2_compress_cmd = "%s c {params} {infile} {outfile}" % lib2_exec
lib2_decompress_cmd = "%s d {infile} {outfile}" % lib2_exec

lib3_compress_cmd = "%s c {params} {infile} {outfile}" % lib3_exec
lib3_decompress_cmd = "%s d {infile} {outfile}" % lib3_exec

diff_cmd = "diff -q {infile} {outfile}"


class Compressor:
    LIB1 = "LIB1"
    LIB2 = "LIB2"
    LIB3 = "LIB3"
    DSRC = "DSRC"

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
        raise RunException("Error running command: %s (exit status: %d)" % (cmd_, ret))

def perform_test(infile_, compressor_, decompressor_, lossy_ = False):
    dsrc_tmp_file = "__out.dsrc"
    fastq_tmp_file = "__out.fastq"
    
    # set-up : cleanup
    if os.path.isfile(dsrc_tmp_file):
        os.remove(dsrc_tmp_file)
        
    if os.path.isfile(fastq_tmp_file):
        os.remove(fastq_tmp_file)

    # this can be done with dictionary...
    if compressor_ == Compressor.LIB1:
        comp_cmd = lib1_compress_cmd
    elif compressor_ == Compressor.LIB2:
        comp_cmd = lib2_compress_cmd
    elif compressor_ == Compressor.LIB3:
        comp_cmd = lib3_compress_cmd
    else:
        comp_cmd = dsrc_compress_cmd

    if decompressor_ == Compressor.LIB1:
        decomp_cmd = lib1_decompress_cmd
    elif decompressor_ == Compressor.LIB2:
        decomp_cmd = lib2_decompress_cmd
    elif decompressor_ == Compressor.LIB3:
        decomp_cmd = lib3_decompress_cmd
    else:
        decomp_cmd = dsrc_decompress_cmd

    if lossy_:
        if compressor_ == Compressor.DSRC:
            params = "-m2 -l -f1,2"
        else:
            params = "L"
    else:
        if compressor_ == Compressor.DSRC:
            params = "-m0"
        else:
            params = ""

    try:
        print "Compressing..."
        run_cmd(comp_cmd.format(params=params,
                                infile=infile_,
                                outfile=dsrc_tmp_file))

        print "Decompressing..."
        run_cmd(decomp_cmd.format(infile=dsrc_tmp_file,
                                  outfile=fastq_tmp_file))  

        if not lossy_:
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
    
    compressors = [Compressor.LIB1, Compressor.LIB2, Compressor.LIB3, Compressor.DSRC]

    # simple test : only lossless mode
    #
    for ci in compressors:
        for co in compressors:
            print "**** Running case: %s->%s ****" % (ci, co)
            perform_test(infile_, ci, co)
    

    # TODO: lossy mode
    #

    
if __name__ == "__main__":
    if len(sys.argv) != 2 or "-h" in sys.argv[1]:
        print "usage: test_dsrc_bin.py <FASTQ filename>"
        exit(0)
        
    run_tests(sys.argv[1])
