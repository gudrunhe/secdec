
CPP_TESTS := $(patsubst %.cpp,%,$(wildcard *.cpp))
CU_TESTS := $(patsubst %.cu,%,$(wildcard *.cu))

# get path to the top level directory
TOPDIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# python executable
PYTHON ?= python3

# c++ compiler
CXX ?= g++

# c++ compiler flags
SECDEC_CONTRIB = $(shell ${PYTHON} -m pySecDecContrib --dirname)
CXXFLAGS = -std=c++17 -O2 -I. -I$(TOPDIR) -I$(SECDEC_CONTRIB)/include

# linker flags
LDFLAGS = -L$(SECDEC_CONTRIB)/lib
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
	rm -rf *.tmp *.log $(CPP_TESTS) $(CU_TESTS)
