DIST_FILENAME=Hydro
DIST_TYPE=app

include ../GLApp/Makefile.mk

# me trying to figure out if numerical inaccuracy compared to javascript version (which uses doubles) is due to build options
#CFLAGS+= -I../GLApp/include -I../TensorMath/include -m32 -mfpmath=387 -mno-sse
#LDFLAGS+= -m32
# default:
INCLUDE+=../GLApp/include
INCLUDE+=../TensorMath/include
INCLUDE+=../Profiler/include
LIBPATHS+=../Profiler/dist/$(PLATFORM)/$(BUILD)
LIBS+=Profiler

$(DIST):: $(OBJPATHS)
	install_name_tool -change dist/$(PLATFORM)/$(BUILD)/libGLApp.dylib ../GLApp/dist/$(PLATFORM)/$(BUILD)/libGLApp.dylib $@
	install_name_tool -change dist/$(PLATFORM)/$(BUILD)/libProfiler.dylib ../Profiler/dist/$(PLATFORM)/$(BUILD)/libProfiler.dylib $@

