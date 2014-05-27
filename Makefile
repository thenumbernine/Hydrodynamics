DIST_FILENAME=Hydro
DIST_TYPE=app

include ../GLApp/Makefile.mk

# me trying to figure out if numerical inaccuracy compared to javascript version (which uses doubles) is due to build options
#CFLAGS_BASE+= -I../GLApp/include -I../TensorMath/include -m32 -mfpmath=387 -mno-sse
#LDFLAGS_BASE+= -m32
# default:
CFLAGS_BASE+= -I../GLApp/include
CFLAGS_BASE+= -I../TensorMath/include
CFLAGS_BASE+= -I../Profiler/include
LDFLAGS_BASE+= -L../Profiler/dist/$(PLATFORM)/$(BUILD) -lProfiler

$(DIST):: $(OBJPATHS)
	install_name_tool -change dist/$(PLATFORM)/$(BUILD)/libGLApp.dylib ../GLApp/dist/$(PLATFORM)/$(BUILD)/libGLApp.dylib $@
	install_name_tool -change dist/$(PLATFORM)/$(BUILD)/libProfiler.dylib ../Profiler/dist/$(PLATFORM)/$(BUILD)/libProfiler.dylib $@

