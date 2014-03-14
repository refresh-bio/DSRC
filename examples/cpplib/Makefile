all: example1 example2

INC_PATH = ../../include/dsrc
LIB_PATH = ../../lib/
DSRC_LIB = -ldsrc

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I$(INC_PATH)

example1: example1.o
	$(CXX) $(CXXFLAGS) -o $@ $? -L$(LIB_PATH) $(DSRC_LIB) $(DEP_LIBS)
	strip $@

example2: example2.o
	$(CXX) $(CXXFLAGS) -o $@ $? -L$(LIB_PATH) $(DSRC_LIB) $(DEP_LIBS)
	strip $@

clean:
	-rm *.o
	-rm example1
	-rm example2
