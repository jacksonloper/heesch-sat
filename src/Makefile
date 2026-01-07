# Check for platform
UNAME_S := $(shell uname -s)

# Build mode: set DEBUG=1 for debug build with symbols and assertions
# Example: make DEBUG=1 render_witness
ifdef DEBUG
    OPT = -g -O0
else
    OPT = -O3 -DNDEBUG
endif

ifeq ($(UNAME_S),Darwin)
    # macOS settings
    INCLUDES = -I/usr/local/include/cryptominisat5 -I/usr/local/include -I/opt/local/include -I/opt/local/include/cairo
    CXX = clang++
    CXXFLAGS = -std=c++17 -Wall -MMD $(INCLUDES) $(OPT) -stdlib=libc++ \
        -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
    LIBS = -L/usr/local/lib -rpath /usr/local/lib -lcryptominisat5
    CAIROLIBS = -L/opt/local/lib -lcairo
else
    # Linux settings
    # Include both /usr/local (for custom LARGEMEM build) and /usr (for system package)
    INCLUDES = -I/usr/local/include/cryptominisat5 -I/usr/include/cryptominisat5 $(shell pkg-config --cflags cairo 2>/dev/null)
    CXX = g++
    CXXFLAGS = -std=c++17 -Wall -MMD $(INCLUDES) $(OPT)
    # Link against /usr/local first (for custom LARGEMEM build), then fall back to system /usr/lib
    LIBS = -L/usr/local/lib -L/usr/lib -Wl,-rpath,/usr/local/lib -lcryptominisat5
    CAIROLIBS = $(shell pkg-config --libs cairo 2>/dev/null)
endif

OBJECTS = sat.o viz.o surrounds.o gen.o report.o render_witness.o
DEPENDS = ${OBJECTS:.o=.d}

all: sat viz gen surrounds report render_witness

sat: sat.o
	$(CXX) -o sat sat.o $(LIBS)

viz: viz.o
	$(CXX) -o viz viz.o $(LIBS) $(CAIROLIBS)

gen: gen.o
	$(CXX) -o gen gen.o $(LIBS)

surrounds: surrounds.o
	$(CXX) -o surrounds surrounds.o $(LIBS)

report: report.o
	$(CXX) -o report report.o $(LIBS)

render_witness: render_witness.o
	$(CXX) -o render_witness render_witness.o $(LIBS)

# Debug version of render_witness with debug symbols for stack traces
# Usage: make render_witness_debug
render_witness_debug: render_witness.cpp
	$(CXX) -std=c++17 -Wall -MMD $(INCLUDES) -g -O0 -o render_witness_debug render_witness.cpp $(LIBS)

# tile: tile.o
#	$(CXX) $(LIBS) -o tile tile.o

-include ${DEPENDS}

.PHONY: clean

clean:
	rm -f ${OBJECTS} ${DEPENDS} sat viz gen surrounds report render_witness render_witness_debug
