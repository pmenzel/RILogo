CXX = g++
CXXFLAGS = -Wall -ansi -pedantic -O3 -g
LDFLAGS = -Wl,-rpath,/usr/local/lib
OBJECTS = ConfigFile.o rilogo.o

all: RILogo makefile

RILogo: makefile $(OBJECTS) 
	$(CXX) -o RILogo $(LDFLAGS) $(OBJECTS)

rilogo.o: rilogo.cpp
	$(CXX) $(CXXFLAGS) -c rilogo.cpp

ConfigFile.o: ConfigFile.cpp ConfigFile.h
	$(CXX) $(CXXFLAGS) -c ConfigFile.cpp

clean: 
	rm -f -v *.o RILogo

.PHONY: clean
