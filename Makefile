#
# Makefile for SOL_Simulations

# ----- Make Macros -----

CXX             = g++
CXXFLAGS        = -g -pedantic -w -Wall -Wextra -std=c++11

TARGETS = main
OBJECTS = Orbit.o Vector.o Particle.o main.o
SRCDIR  = src/

# ----- Make rules -----

all:	$(TARGETS)

clean:
	rm -rf $(TARGETS) $(OBJECTS)
	
main:	$(OBJECTS)
	$(CXX) -o main $(OBJECTS)

Orbit.o:	$(SRCDIR)Orbit.cc $(SRCDIR)Orbit.h $(SRCDIR)Vector.h $(SRCDIR)Particle.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Orbit.cc

Vector.o:	$(SRCDIR)Vector.h $(SRCDIR)Vector.cc
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Vector.cc

Particle.o: $(SRCDIR)Particle.h $(SRCDIR)Particle.cc $(SRCDIR)Vector.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)Particle.cc

main.o: $(SRCDIR)main.cc $(SRCDIR)Orbit.h $(SRCDIR)Vector.h $(SRCDIR)Particle.h
	$(CXX) -c $(CXXFLAGS) $(SRCDIR)main.cc
