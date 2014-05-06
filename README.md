# DSRC

  DSRC is a toolkit designed for efficient high-performance compression of sequencing reads stored in FASTQ format, where it's main features are:
* Effective multithreaded compression of FASTQ files.
* Full support for Illumina, ABI SOLiD, and 454/Ion Torrent dataset formats with non-standard (AGCTN) IUPAC base values.
* Support for lossy quality values compression using Illumina binning scheme.
* Support for lossy IDs compression keeping only key fields selected by user.
* Pipes support for easy integration with current pipelines.
* Python and C++ libraries allowing to integrate DSRC archives in own applications.
* Availability for Linux, Mac OSX and Windows 64-bit operating systems.
* Open source C++ code under GNU GPL 2 license.

For more information please check out the [official website](http://sun.aei.polsl.pl/dsrc/).


# Building

## Pre-built binaries and libraries

Pre-built binaries and libraries for Linux, Mac OSX and Windows platforms can be downloaded from the [official website](http://sun.aei.polsl.pl/dsrc/download.html).


## Build prerequisites

### Linux

DSRC binaries and C++ library can be compiled in two ways, depending on the selection of multithreading support library - for each a different makefile file is provided. In the first case, _boost::threads_ library will be used, which is needed to be present on the build system. In the second - _g++_ compiler with c++11 support (version >= 4.8).

By default, binaries and libraries are compiled using _g++_, however compiling using _Clang_ or _Intel icpc_ should also succeeed without any problems.


### Mac OSX

On Mac OSX _Clang_ compiler will be used with c++11 support, so make sure to have _Clang_ in version >= 3.3 installed.


### Windows

To compile DSRC under Windows OS, _Microsoft Visual Studio_ 2010 or 2012 is required. DSRC binaries and C++ library can be compiled in two ways, depending on the selection of multithreading support library - for each a different VS solution file is provided. When compiling using VS2010 the _boost::threads_ library will be used to provide multithreading support, so make sure to have _boost::threads_ library installed and _boost_ library paths properly configured in Visual Studio. In case of using VS2012 c++11 standard implementation will be used to provide threading support.

There should be also no problems when compiling DSRC using _MinGW-32-x64_ with provided Makefile files.


### Python library

To build DSRC Python library, _boost::python_ library in development version and _boost::build_ tool _bjam_ are need to be present on the system. Next, in the _Jamroot_ configuration file in _py_ directory a local _boost_ installation directory needs to be specified:

```
# To compile DSRC Python module please specify your boost installation directory below
#
use-project boost 
	: /absolute/path/to/boost/directory/ ;
```

Python library will be built using a default compilation toolset available on the build platform (auto selected by _bjam_), however in order to specify a different one append

    <toolset>name

to the compilation flags as exmplained in the _Jamroot_ file

```
# Specify toolset according to your platform manually in case of compilation problems in form: '<toolset>gcc'
# Available toolsets:
#	- Windows: msvc-*
#	- Linux: gcc, clang
#	- Mac OSX: darwin, gcc
	: <variant>release <address-model>64 <link>shared <runtime-link>shared <debug-symbols>off <inlining>full <optimization>speed <warnings>on <cxxflags>"-O2 -m64 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DUSE_BOOST_THREAD" ;
```

## Building on Linux

### Binary

To compile DSRC using _boost::threads_ with static linking, in the main directory type:

    make bin

To compile DSRC using _g++_ >= 4.8 with c++11 standard and dynamic linking:
    
    make -f Makefile.c++11 bin

The resulting _dsrc_ binary will be placed in _bin_ subdirectory.



### C++ library

To compile C++ DSRC library using _boost::threads_:

    make lib

To compile DSRC using _g++_ >= 4.8 with c++11:
    
    make -f Makefile.c++11 lib

The resulting _libdsrc.a_ library will be placed in _lib_ subdirectory.


### Python library

To compile DSRC Python library:

    make pylib
    
The resulting _pydsrc.so_ library will be available in _py_ subdirectory.



## Building on Mac OSX

### Binary

To compile DSRC binary, in the main directory type:

    make -f Makefile.osx bin
        
The resulting _dsrc_ binary will be placed in _bin_ subdirectory.


### C++ library

To compile DSRC C++ library:
    
    make -f Makefile.osx lib

The resulting _libdsrc.a_ library will be placed in _lib_ subdirectory.


### Python library

To compile DSRC Python library:

    make -f Makefile.osx pylib
    
The resulting _pydsrc.so_ library will be available in _py_ subdirectory.



## Building on Windows 64-bit

### Binary

To compile DSRC using _Visual Studio 2010_ with _boost::threads_ for multithreading support use the _dsrc20-vs2k10.sln_ solution file. However, to compile DSRC using _Visual Studio 2012_ with c++11 threads use the _dsrc20-vs2k12.sln_. 

To compile DSRC executable, select `Release|x64` configuration and build.

The resulting _dsrc.exe_ executable will be placed in _bin_ subdirectory.


### C++ library

To compile DSRC library, select `Release Lib|x64` configuration and build.

The resulting _dsrc.lib_ library will be placed in _lib_ subdirectory.


### Python library

To compile DSRC Python library in the _py_ subdirectory type:

    bjam
    
The resulting _pydsrc.pyd_ library will be available in _py_ subdirectory.



# Usage

DSRC can be run from the command prompt:

    dsrc <c|d> [options] <input_file_name> <output_file_name>

in one of two modes:
* `c` — compression,
* `d` — decompression.

## Available options

### Compression options
* `-d<n>` — DNA compression mode: `0–3`, default: `0`
* `-q<n>` — Quality compression mode: `0–2`, default: `0`
* `-f<1,...>` — keep only those fields no. in ID field string, default: ` ` (keep all)
* `-b<n>` — FASTQ input buffer size in MB, default: `8`
* `-m<n>` — Automated compression mode (one of the three preset combination of other pa-
rameters): `0–2`
* `-o<n>` — Quality offset, 0 for auto selection, default: `0`
* `-l` — use Quality lossy mode (Illumina binning scheme), default: `false`
* `-c` — calculate and check CRC32 checksum calculation per block (slows the compression
about twice), default: `false`

### Automated compression modes
* `-m0` — fast mode, equivalent to: `-d0 -q0 -b8`
* `-m1` — medium mode, equivalent to: `-d2 -q2 -b64`
* `-m2` — best mode, equivalent to: `-d3 -q2 -b256`

### Options for both compression and decompression
* `-t<n>` — processing threads number, default: max available hardware threads
* `-s` — use stdin/stdout for reading/writing raw FASTQ files data (stderr is used for info/warning
messages)


## Usage examples
Compress `SRR001471.fastq` file saving DSRC archive to `SRR001471.dsrc`:

    dsrc c SRR001471.fastq SRR001471.dsrc
    
Compress file in the fast mode with CRC32 checking and using `4` threads:

    dsrc c -m0 -c -t4 SRR001471.fastq SRR001471.dsrc
    
Compress file using DNA and Quality compression level `2` and using `512` MB buffer:

    dsrc c -d2 -q2 -b512 SRR001471.fastq SRR001471.dsrc
    
Compress file in the best mode with lossy Quality mode and preserving only `1–4` fields from
record IDs:

    dsrc c -m2 -l -f1,2,3,4 SRR001471.fastq SRR001471.dsrc
    
Compress in the best mode reading raw FASTQ data from stdin:

    cat SRR001471.fastq | dsrc c -m2 -s SRR001471.dsrc
    
Decompress `SRR001471.dsrc` archive saving output FASTQ file to `SRR001471.out.fastq`:

    dsrc d SRR001471.dsrc SRR001471.out.fastq

Decompress archive using `4` threads and streaming raw FASTQ data to stdout:

    dsrc d -t4 -s SRR001471.dsrc > SRR001471.out.fastq
    
