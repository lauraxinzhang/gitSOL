#
# Makefile for SOL_Simulations

# ----- Make Macros -----

CXX             = g++
CXXFLAGS        = -g -pedantic -w -fopenmp -Wall -Wextra -std=c++11 -O3
TARGETS = main
OBJECTS = Orbit.o Orbit.test.o Orbit.write.o Mirror.o Vector.o Particle.o Pusher.o main.o
SRCDIR  = src/

# ----- Make rules -----

all:	$(TARGETS)

clean:
	rm -rf $(TARGETS) $(OBJECTS)
	
main:	$(OBJECTS)
	$(CXX) $(CXXFLAGS) -o main $(OBJECTS) $(LDFLAGS)

Orbit.o:	$(SRCDIR)Orbit.cc $(SRCDIR)Orbit.h $(SRCDIR)Vector.h $(SRCDIR)Particle.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Orbit.cc $(LDFLAGS)

Orbit.test.o:	$(SRCDIR)Orbit.test.cc $(SRCDIR)Orbit.h $(SRCDIR)Vector.h $(SRCDIR)Particle.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Orbit.test.cc

Orbit.write.o:	$(SRCDIR)Orbit.write.cc $(SRCDIR)Orbit.h $(SRCDIR)Vector.h $(SRCDIR)Particle.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Orbit.write.cc

Mirror.o:	$(SRCDIR)Mirror.cc $(SRCDIR)Mirror.h $(SRCDIR)Vector.h $(SRCDIR)Particle.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Mirror.cc $(LDFLAGS)

Vector.o:	$(SRCDIR)Vector.h $(SRCDIR)Vector.cc
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Vector.cc

Particle.o: $(SRCDIR)Particle.h $(SRCDIR)Particle.cc $(SRCDIR)Vector.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Particle.cc

main.o: $(SRCDIR)main.cc $(SRCDIR)Pusher.h $(SRCDIR)Orbit.h $(SRCDIR)Mirror.h $(SRCDIR)Vector.h $(SRCDIR)Particle.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)main.cc $(LDFLAGS)
