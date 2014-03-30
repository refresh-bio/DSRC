# make sure to configure properly Jamroot file and have boost libraries installed
.PHONY: all clean

PYLIB_NAME = pydsrc.so
BJAM_BUILD_ROOT = bin
BJAM_BUILD_PATH = $(BJAM_BUILD_ROOT)/*/release/address-model-64

all:
	bjam
	mv $(BJAM_BUILD_PATH)/$(PYLIB_NAME) ./

clean:
	-rm $(PYLIB_NAME)
	-rm -r $(BJAM_BUILD_ROOT)
