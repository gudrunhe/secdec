
TESTS := $(patsubst %.cpp,%,$(wildcard *.cpp))

# get path to the top level directory
TOPDIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# python executable
PYTHON ?= python

# c++ compiler
CXX ?= g++

# c++ compiler flags
CXXFLAGS = -std=c++11 -O2 -I. -I$(TOPDIR)
ifneq (x$(SECDEC_CONTRIB), x)
CXXFLAGS += -I$(SECDEC_CONTRIB)/include
endif

# linker flags
ifeq (x$(SECDEC_CONTRIB), x)
LDFLAGS =
else
LDFLAGS = -L$(SECDEC_CONTRIB)/lib
endif
LDFLAGS += -L$(TOPDIR) -lcuba -lm

.PHONY: clean doctest $(TESTS)

doctest: $(TESTS)

$(TESTS):
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -fPIC $@.cpp -o $@
	./$@ > $@.tmp
	diff $@.tmp $@.txt | tee $@.log
	test ! -s $@.log 

clean: 
	rm -rf *.tmp *.log $(TESTS)