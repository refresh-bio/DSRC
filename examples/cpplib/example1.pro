TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    example1.cpp

QMAKE_CXXFLAGS += -std=c++11

INCLUDEPATH += ../../include/dsrc

# Specify path to DSRC library file
LIBS += $$_PRO_FILE_PWD_/../../lib/libdsrc.a

LIBS += -lpthread
#LIBS += -lboost_thread
