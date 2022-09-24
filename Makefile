DIST_FILENAME=Hydro
DIST_TYPE=app

include ../Common/Base.mk
include ../GLApp/Include.mk
include ../Tensor/Include.mk
include ../Profiler/Include.mk
include ../Parallel/Include.mk

# me trying to figure out if numerical inaccuracy compared to javascript version (which uses doubles) is due to build options
#CXXFLAGS+= -m32 -mfpmath=387 -mno-sse
#LDFLAGS+= -m32
