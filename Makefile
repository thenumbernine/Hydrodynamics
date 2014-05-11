DISTDIR=dist/osx
DIST=$(DISTDIR)/hydro

SOURCES=$(shell find src -type f)
HEADERS=$(shell find include -type f)
OBJECTS=HydroApp.o

OBJDIR=obj/osx/release
OBJPATHS=$(addprefix $(OBJDIR)/, $(OBJECTS))

CC=clang++
CFLAGS_BASE=-c -Wall -std=c++11 -Iinclude -I../GLApp/include -I../TensorMath/include
CFLAGS_DEBUG=$(CFLAGS_BASE) -O0 -mfix-and-continue -gdwarf-2
CFLAGS_RELEASE=$(CFLAGS_BASE) -Os
CFLAGS=$(CFLAGS_RELEASE)

LD=clang++
LDFLAGS_BASE=-L../GLApp/dist/osx -lGLApp -lSDL -lSDLmain -framework Cocoa -framework OpenGL
LDFLAGS_DEBUG=$(LDFLAGS_BASE)
LDFLAGS_RELEASE=$(LDFLAGS_BASE)
LDFLAGS=$(LDFLAGS_RELEASE)

.PHONY: all prep clean distclean test

all: prep $(DIST)

prep:
	-mkdir -p $(OBJDIR)
	-mkdir -p $(DISTDIR)

$(OBJDIR)/%.o : src/%.cpp $(HEADERS)
	-mkdir -p $(OBJDIR)/$(<D)
	$(CC) $(CFLAGS) -o $@ $<

$(DIST): $(OBJPATHS)
	$(LD) $(LDFLAGS) -o $@ $^
	install_name_tool -change dist/osx/libGLApp.dylib ../GLApp/dist/osx/libGLApp.dylib dist/osx/hydro

clean:
	-rm -f $(OBJPATHS)

distclean:
	-rm -f $(DIST)

test: $(DIST)
	$(MAKE) -C test run

