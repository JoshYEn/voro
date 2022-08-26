# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Load the common configuration file
include ./config.mk

# List of executables
EXECUTABLES=test_htr

# Makefile rules
all: $(EXECUTABLES)

test_htr: htr.cc
	$(CXX) $(CFLAGS) $(T_INC) $(T_LIB) -o test_htr.o htr.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
