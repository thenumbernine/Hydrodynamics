# PLATFORM: osx
# BUILD: debug, release

OBJECTS=HydroApp.o Profile.o
SOURCES=$(shell find src -type f)
HEADERS=$(shell find include -type f)

DISTDIR_BASE=dist
DIST_FILENAME=hydro
DISTDIR=$(DISTDIR_BASE)/$(PLATFORM)/$(BUILD)
DIST=$(DISTDIR)/$(DIST_FILENAME)

OBJDIR_BASE=obj
OBJDIR=$(OBJDIR_BASE)/$(PLATFORM)/$(BUILD)
OBJPATHS=$(addprefix $(OBJDIR)/, $(OBJECTS))

CC=clang++
CFLAGS_BASE=-c -Wall -std=c++11 -Iinclude -I../GLApp/include -I../TensorMath/include
CFLAGS_debug=-O0 -mfix-and-continue -gdwarf-2 -DDEBUG
CFLAGS_release=-O3 -DNDEBUG

LD=clang++
LDFLAGS_BASE=-L../GLApp/dist/$(PLATFORM)/$(BUILD) -lGLApp -lSDL -lSDLmain -framework Cocoa -framework OpenGL

.PHONY: default
default: osx

.PHONY: help
help:
	echo "make <platform>"
	echo "platform: osx"

.PHONY: osx
osx:
	$(MAKE) PLATFORM="osx" build_platform

.PHONY: build_platform
build_platform: $(PLATFORM)_debug $(PLATFORM)_release

.PHONY: $(PLATFORM)_debug
$(PLATFORM)_debug:
	$(MAKE) BUILD="debug" dist

.PHONY: $(PLATFORM)_release
$(PLATFORM)_release:
	$(MAKE) BUILD="release" dist

.PHONY: dist
dist: CFLAGS= $(CFLAGS_BASE)
dist: CFLAGS+= $(CFLAGS_$(BUILD))
dist: LDFLAGS= $(LDFLAGS_BASE)
dist: LDFLAGS+= $(LDFLAGS_$(BUILD))
dist: $(DIST)

$(OBJDIR)/%.o : src/%.cpp $(HEADERS)
	-mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $<

$(DIST): $(OBJPATHS)
	-mkdir -p $(@D)
	$(LD) $(LDFLAGS) -o $@ $^
	install_name_tool -change dist/$(PLATFORM)/$(BUILD)/libGLApp.dylib ../GLApp/dist/$(PLATFORM)/$(BUILD)/libGLApp.dylib $@

.PHONY: clean
clean:
	-rm -fr $(OBJDIR_BASE)

.PHONY: distclean
distclean:
	-rm -fr $(DISTDIR_BASE)

.PHONY: test
test: $(DIST)
	$(MAKE) -C test run

