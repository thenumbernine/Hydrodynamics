# PLATFORM: osx
# BUILD: debug, release

DIST_FILENAME=hydro
DIST_TYPE=app

OBJECTS=HydroApp.o Profile.o
SOURCES=$(shell find src -type f)
HEADERS=$(shell find include -type f)

include ../GLApp/Makefile.mk

CFLAGS_BASE+=-I../GLApp/include -I../TensorMath/include

$(DIST):: $(OBJPATHS)
	install_name_tool -change dist/$(PLATFORM)/$(BUILD)/libGLApp.dylib ../GLApp/dist/$(PLATFORM)/$(BUILD)/libGLApp.dylib $@

