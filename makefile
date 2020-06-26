IDIR=${CURDIR}/src

	
program=Fourier3D

SOURCES := $(shell find $(IDIR) -name '*.cpp')

CXX = g++ -O3 -s -DNDEBUG
#CXX = g++ -g -W -Wall -Werror
#CXX = g++ -O0 -g3 -Wall -c
CXXFLAGS = -std=c++0x
LDFLAGS = -lfftw3 -lfftw3f -L./lib/ -L./
#LDFLAGS += -g

OBJECTS = $(SOURCES:.cpp=.o)

all: ${program}

build: ${program}

debug: CXX = g++ -g -W -Wall -Werror 
debug: LDFLAGS += -g
debug: ${program}

 
${program}: CXXFLAGS += $(foreach d, $(includepath), -I$d)
${program}: LDFLAGS += $(foreach d, $(libpath), -L$d)
${program}: $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^

clean:
	rm -f src/*.o ${program} src/*.d

.cpp.o:
	$(CXX) -MD -MP $(CXXFLAGS) -o $@ -c $< $(LDFLAGS)
