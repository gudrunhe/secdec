
CPP_TESTS := $(patsubst %.cpp,%,$(wildcard *.cpp))
CU_TESTS := $(patsubst %.cu,%,$(wildcard *.cu))

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
ifneq (x$(CXX), xnvcc)
CXXFLAGS += -x c++
LDFLAGS += -pthread
endif
LDFLAGS += -L$(TOPDIR) -lcuba -lgsl -lgslcblas -lm

.PHONY: clean doctest $(CPP_TESTS) $(CU_TESTS)

doctest: $(CPP_TESTS) $(CU_TESTS)

$(CPP_TESTS):
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)
	./$@ > $@.tmp
	diff $@.tmp $@.txt | tee $@.log
	test ! -s $@.log

$(CU_TESTS):
	$(CXX) $(CXXFLAGS) $@.cu -o $@ $(LDFLAGS)
	./$@ > $@.tmp
	diff $@.tmp $@.txt | tee $@.log
	test ! -s $@.log

clean:
	rm -rf *.tmp *.log $(TESTS)
