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
    INCLUDES = -I/usr/include/cryptominisat5 $(shell pkg-config --cflags cairo)
    CXX = g++
    OPT = -O3 -DNDEBUG
    CXXFLAGS = -std=c++17 -Wall -MMD $(INCLUDES) $(OPT)
    LIBS = -lcryptominisat5
    CAIROLIBS = $(shell pkg-config --libs cairo)
endif

OBJECTS = sat.o viz.o surrounds.o gen.o report.o render_witness.o
DEPENDS = ${OBJECTS:.o=.d}

all: sat viz gen surrounds report render_witness

sat: sat.o 
	$(CXX) $(LIBS) -o sat sat.o

viz: viz.o
	$(CXX) $(LIBS) $(CAIROLIBS) -o viz viz.o

render_witness: render_witness.o
	$(CXX) -o render_witness render_witness.o $(LIBS) 

gen: gen.o
	$(CXX) $(LIBS) -o gen gen.o

surrounds: surrounds.o 
	$(CXX) $(LIBS) -o surrounds surrounds.o

report: report.o 
	$(CXX) $(LIBS) -o report report.o

# tile: tile.o 
#	$(CXX) $(LIBS) -o tile tile.o

-include ${DEPENDS}

.PHONY: clean

clean:
	rm ${OBJECTS} ${DEPENDS}
