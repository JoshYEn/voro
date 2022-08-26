# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Load the common configuration file
include ./config.mk

# List of executables
EXECUTABLES=test_vf

# Makefile rules
all: $(EXECUTABLES)

test_vf: htr_vf.cc
	$(CXX) $(CFLAGS) $(T_INC) $(T_LIB) -o test_vf.o htr_vf.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
