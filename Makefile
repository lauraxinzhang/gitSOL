#
# Makefile for SOL_Simulations

# ----- Make Macros -----

CXX             = g++
CXXFLAGS        = -g -pedantic -w -fopenmp -Wall -Wextra -std=c++11 -O3
TARGETS = main
OBJECTS = Orbit.o Orbit.test.o Orbit.write.o Orbit.thermal.o Struct.o Mirror.o Vector.o Particle.o main.o
SRCDIR  = src/

# ----- Make rules -----

all:	$(TARGETS)

clean:
	rm -rf $(TARGETS) $(OBJECTS)
	
main:	$(OBJECTS)
	$(CXX) $(CXXFLAGS) -o main $(OBJECTS) $(LDFLAGS)

Orbit.o:	$(SRCDIR)Orbit.cc $(SRCDIR)Orbit.h $(SRCDIR)Vector.h $(SRCDIR)Particle.h $(SRCDIR)Struct.h $(SRCDIR)Constants.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Orbit.cc $(LDFLAGS)

Orbit.test.o:	$(SRCDIR)Orbit.test.cc $(SRCDIR)Orbit.h $(SRCDIR)Vector.h $(SRCDIR)Particle.h  $(SRCDIR)Constants.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Orbit.test.cc

Orbit.write.o:	$(SRCDIR)Orbit.write.cc $(SRCDIR)Orbit.h $(SRCDIR)Vector.h $(SRCDIR)Particle.h $(SRCDIR)Constants.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Orbit.write.cc

Orbit.thermal.o:$(SRCDIR)Orbit.thermal.cc $(SRCDIR)Orbit.h $(SRCDIR)Vector.h $(SRCDIR)Particle.h $(SRCDIR)Constants.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Orbit.thermal.cc

Struct.o: $(SRCDIR)Struct.cc $(SRCDIR)Orbit.h $(SRCDIR)Struct.h $(SRCDIR)Vector.h $(SRCDIR)Particle.h $(SRCDIR)Constants.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Struct.cc

Mirror.o:	$(SRCDIR)Mirror.cc $(SRCDIR)Mirror.h $(SRCDIR)Vector.h $(SRCDIR)Particle.h $(SRCDIR)Constants.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Mirror.cc $(LDFLAGS)

Vector.o:	$(SRCDIR)Vector.h $(SRCDIR)Vector.cc
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Vector.cc

Particle.o: $(SRCDIR)Particle.h $(SRCDIR)Particle.cc $(SRCDIR)Vector.h $(SRCDIR)Constants.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Particle.cc

main.o: $(SRCDIR)main.cc $(SRCDIR)Pusher.h $(SRCDIR)Orbit.h $(SRCDIR)Mirror.h $(SRCDIR)Vector.h $(SRCDIR)Particle.h $(SRCDIR)Struct.h  $(SRCDIR)Constants.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)main.cc $(LDFLAGS)
