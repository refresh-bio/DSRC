all: lib bin clean

APP_NAME = dsrc
LIB_NAME = libdsrc.a

INC_DIR = ../include/dsrc
AR_BIN = ar
AR_FLAGS = -cvq

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

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

bin: main.o $(APP_OBJS)
	$(CXX) $(CXXFLAGS) -o $(APP_NAME) main.o $(APP_OBJS) $(DEP_LIB_DIRS) $(DEP_LIBS)
	strip $(APP_NAME)

lib: $(APP_OBJS) $(LIB_OBJS)
	ar $(AR_FLAGS) $(LIB_NAME) $(APP_OBJS) $(LIB_OBJS)

clean:
	-rm *.o
	-rm $(APP_NAME)
	-rm $(LIB_NAME)
