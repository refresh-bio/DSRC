.PHONY: all lib bin examples clean

all: lib bin examples

export

CXX = g++
CXXFLAGS = -O2 -m64
CXXFLAGS += -DNDEBUG -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE
CXXFLAGS += -Wall

# there are problems when statically linking with pthreads and libstdc++
# so use only shared linking
#CXXFLAGS += -static

# compile using c++11
CXXFLAGS += -std=c++11

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

clean:
	cd src; ${MAKE} clean
	cd examples/cpplib; ${MAKE} clean
	-rm -r $(LIB_DIR)
	-rm -r $(BIN_DIR)