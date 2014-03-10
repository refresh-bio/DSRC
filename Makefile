.PHONY: all bin lib pylib examples clean

all: bin lib examples

export

CXX = g++
CXXFLAGS = -O2 -m64
CXXFLAGS += -DNDEBUG -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE
CXXFLAGS += -Wall #-pedantic

# comment the line below to compile and link in shared mode
# CXXFLAGS += -static

# compile using boost::thread 
# boost::thread from 1.50+ explicitely requires boost::system library
CXXFLAGS += -DUSE_BOOST_THREAD
DEP_LIBS += -lboost_thread -lboost_system

# necessary library to link
# (even when using boost, remember to link it _after_ linking with boost::thread)
DEP_LIBS += -lpthread


bin:
	cd src; ${MAKE} bin

lib:
	cd src; ${MAKE} lib

pylib:
	# make sure to configure properly Jamroot file and have boost libraries installed
	cd py; bjam

examples:
	cd examples/cpplib; ${MAKE}

clean:
	cd src; ${MAKE} clean
	cd examples/cpplib; ${MAKE} clean
