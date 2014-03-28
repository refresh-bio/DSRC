TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    example2.cpp

INCLUDEPATH += ../../include/dsrc

# Specify path to DSRC library file
LIBS += /home/lucas/dev/workspace/dsrc20-dev/lib/libdsrc.a

LIBS += -lpthread
LIBS += -lboost_thread
