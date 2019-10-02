OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR) -I$(HEEDDIR)

# Debug flags
#CFLAGS += -g

LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += `root-config --glibs` -lGeom -lgfortran -lm
#LDFLAGS += -g

primaries: primaries.C 
	$(CXX) $(CFLAGS) primaries.C
	$(CXX) -o primaries primaries.o $(LDFLAGS)
	rm primaries.o

primaries_single: primaries_single.C 
	$(CXX) $(CFLAGS) primaries_single.C
	$(CXX) -o primaries_single primaries_single.o $(LDFLAGS)
	rm primaries_single.o

