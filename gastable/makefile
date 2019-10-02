INCDIR = $(GARFIELD_HOME)/Include
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CXX = `root-config --cxx`
CFLAGS = `root-config --cflags` \
	-O3 -W -Wall -Wextra -Wno-long-long \
	-fno-common \
	-I$(INCDIR)

LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += `root-config --glibs` -lGeom -lgfortran -lm

gastable: gastable.cc
	$(CXX) $(CFLAGS) -c gastable.cc
	$(CXX) $(CFLAGS) -o gastable gastable.o $(LDFLAGS)
	rm gastable.o

gastable_gen: gastable_gen.cc
	$(CXX) $(CFLAGS) -c gastable_gen.cc
	$(CXX) $(CFLAGS) -o gastable_gen gastable_gen.o $(LDFLAGS)
	rm gastable_gen.o

gastable_read: gastable_read.cc
	$(CXX) $(CFLAGS) -c gastable_read.cc
	$(CXX) $(CFLAGS) -o gastable_read gastable_read.o $(LDFLAGS)
	rm gastable_read.o
