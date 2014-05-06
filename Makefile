DISTDIR=dist/osx
DIST=$(DISTDIR)/hydro

OBJECTS=BoundaryMethod.o EquationOfState.o ExplicitMethod.o FluxMethod.o Hydro.o HydroApp.o InitialConditions.o Solver.o

OBJDIR=obj/osx/release
OBJPATHS=$(addprefix $(OBJDIR)/, $(OBJECTS))

CC=clang++
CFLAGS_BASE=-c -Wall -Iinclude -I../GLApp/include -std=c++11
CFLAGS_DEBUG=$(CFLAGS_BASE) -O0 -mfix-and-continue -gdwarf-2
CFLAGS_RELEASE=$(CFLAGS_BASE) -Os
CFLAGS=$(CFLAGS_DEBUG)

LD=clang++
LDFLAGS_BASE=-L../GLApp/dist/osx -lGLApp -lSDL -lSDLmain -framework Cocoa -framework OpenGL
LDFLAGS_DEBUG=$(LDFLAGS_BASE)
LDFLAGS_RELEASE=$(LDFLAGS_BASE)
LDFLAGS=$(LDFLAGS_DEBUG)

.PHONY: all prep clean distclean test

all: prep $(DIST)

prep:
	-mkdir -p $(OBJDIR)
	-mkdir -p $(DISTDIR)

$(OBJDIR)/%.o : src/%.cpp
	-mkdir -p $(OBJDIR)/$(^D)
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

