.PHONY: all lib bin pylib examples clean

all: bin lib examples

export

CXX = g++
CXXFLAGS = -O2 -m64
CXXFLAGS += -DNDEBUG -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE
CXXFLAGS += -Wall #-pedantic

# comment the line below to compile and link in shared mode
CXXFLAGS += -static

# by default compile using boost::thread 
# boost::thread from 1.50+ explicitely requires boost::system library
CXXFLAGS += -DUSE_BOOST_THREAD
DEP_LIBS += -lboost_thread -lboost_system

# necessary library to link
# (even when using boost, remember to link it _after_ linking with boost::thread)
DEP_LIBS += -lpthread

APP_NAME = dsrc
LIB_NAME = libdsrc.a
LIB_DIR = lib
BIN_DIR = bin
SRC_DIR = src

bin:
	cd src; ${MAKE} bin
	test -d $(BIN_DIR) || mkdir $(BIN_DIR)
	mv $(SRC_DIR)/$(APP_NAME) $(BIN_DIR)/

lib:
	cd src; ${MAKE} lib
	test -d $(LIB_DIR) || mkdir $(LIB_DIR)
	mv $(SRC_DIR)/$(LIB_NAME) $(LIB_DIR)/

examples:
	cd examples/cpplib; ${MAKE}

pylib:
	cd py; ${MAKE}

clean:
	cd src; ${MAKE} clean
	cd examples/cpplib; ${MAKE} clean
	cd py; ${MAKE} clean
	-rm -r $(LIB_DIR)
	-rm -r $(BIN_DIR)
