# Check for platform
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    # macOS settings
    INCLUDES = -I/usr/local/include/cryptominisat5 -I/usr/local/include -I/opt/local/include -I/opt/local/include/cairo
    CXX = clang++
    OPT = -O3 -DNDEBUG
    CXXFLAGS = -std=c++17 -Wall -MMD $(INCLUDES) $(OPT) -stdlib=libc++ \
        -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
    LIBS = -L/usr/local/lib -rpath /usr/local/lib -lcryptominisat5
    CAIROLIBS = -L/opt/local/lib -lcairo
else
    # Linux settings
    INCLUDES = -I/usr/include/cryptominisat5 $(shell pkg-config --cflags cairo 2>/dev/null)
    CXX = g++
    OPT = -O3 -DNDEBUG
    CXXFLAGS = -std=c++17 -Wall -MMD $(INCLUDES) $(OPT)
    LIBS = -lcryptominisat5
    CAIROLIBS = $(shell pkg-config --libs cairo 2>/dev/null)
endif

OBJECTS = sat.o viz.o surrounds.o gen.o report.o
DEPENDS = ${OBJECTS:.o=.d}

all: sat viz gen surrounds report

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

# render_witness is now a Python script (render_witness.py)

# tile: tile.o
#	$(CXX) $(LIBS) -o tile tile.o

-include ${DEPENDS}

.PHONY: clean

clean:
	rm -f ${OBJECTS} ${DEPENDS} sat viz gen surrounds report
