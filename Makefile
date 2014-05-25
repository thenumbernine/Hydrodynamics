DIST_FILENAME=hydro
DIST_TYPE=app

include ../GLApp/Makefile.mk

# me trying to figure out if numerical inaccuracy compared to javascript version (which uses doubles) is due to build options
#CFLAGS_BASE+= -I../GLApp/include -I../TensorMath/include -m32 -mfpmath=387 -mno-sse
#LDFLAGS_BASE+= -m32
# default:
CFLAGS_BASE+= -I../GLApp/include -I../TensorMath/include

$(DIST):: $(OBJPATHS)
	install_name_tool -change dist/$(PLATFORM)/$(BUILD)/libGLApp.dylib ../GLApp/dist/$(PLATFORM)/$(BUILD)/libGLApp.dylib $@

