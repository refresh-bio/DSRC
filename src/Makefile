all: bin lib

INC_DIR = ../include/dsrc
SRC_DIR = .
LIB_DIR = ../lib
BIN_DIR = ../bin

AR_BIN = ar

APP_OBJS = DsrcOperator.o \
	DsrcWorker.o \
	DsrcIo.o \
	DsrcFile.o \
	BlockCompressor.o \
	DnaModelerHuffman.o \
	QualityPositionModeler.o \
	QualityRLEModeler.o \
	TagModeler.o \
	RecordsProcessor.o \
	FastqParser.o \
	FastqIo.o \
	FastqStream.o \
	FileStream.o \
	StdStream.o \
	huffman.o

LIB_OBJS = DsrcArchive.o \
	DsrcModule.o \
	BlockCompressorExt.o \
	FastqFile.o \
	Configurable.o

APP_NAME = dsrc
LIB_NAME = libdsrc.a
AR_FLAGS = -cvq

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

bin: main.o $(APP_OBJS)
	$(CXX) $(CXXFLAGS) -o $(APP_NAME) main.o $(APP_OBJS) $(DEP_LIB_DIRS) $(DEP_LIBS)
	strip $(APP_NAME)
	test -d $(BIN_DIR) || mkdir $(BIN_DIR)
	mv $(APP_NAME) $(BIN_DIR)/


lib: $(APP_OBJS) $(LIB_OBJS)
	ar $(AR_FLAGS) $(LIB_NAME) $(APP_OBJS) $(LIB_OBJS)
	test -d $(LIB_DIR) || mkdir $(LIB_DIR)
	mv $(LIB_NAME) $(LIB_DIR)/

clean:
	-rm *.o
	-rm $(BIN_DIR)/$(APP_NAME)
	-rm $(LIB_DIR)/$(LIB_NAME)
