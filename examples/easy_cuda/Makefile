INTEGRAL = easy

# CXX
export CXX = nvcc

# cuda architecture
export SECDEC_WITH_CUDA = compute_30

# CXXFLAGS
ifeq (x$(SECDEC_CONTRIB), x)
CXXFLAGS =
else
CXXFLAGS = -I$(SECDEC_CONTRIB)/include
endif
CXXFLAGS += -dc -x cu -std=c++11 -O2 -I$(INTEGRAL) -DSECDEC_WITH_CUDA=$(SECDEC_WITH_CUDA)

# LDFLAGS
ifeq (x$(SECDEC_CONTRIB), x)
LDFLAGS =
else
LDFLAGS = -L$(SECDEC_CONTRIB)/lib
endif
LDFLAGS += -L$(INTEGRAL) -l$(INTEGRAL) -lgsl -lgslcblas -lm

# phony targets
.PHONY : all
all : integrate_$(INTEGRAL) pylink

.PHONY : clean
clean :
	rm -f *.o *.so *.a integrate_$(INTEGRAL)
	for dir in *; do if [ -d $$dir ]; then $(MAKE) -C $$dir $@; fi; done

.PHONY : very-clean
very-clean :
	rm -rf */ *.o *.so *.a integrate_$(INTEGRAL)

.PHONY : pylink
pylink : $(INTEGRAL)/lib$(INTEGRAL)_pylink.so

.PHONY : static
static : $(INTEGRAL)/lib$(INTEGRAL).a

integrate_$(INTEGRAL) : integrate_$(INTEGRAL).o $(INTEGRAL)/lib$(INTEGRAL).a
	$(CXX) -o $@ $< $(LDFLAGS)

integrate_$(INTEGRAL).o : integrate_$(INTEGRAL).cpp $(INTEGRAL)/$(INTEGRAL).hpp
	$(CXX) $(CXXFLAGS) -o $@ $<

$(INTEGRAL)/$(INTEGRAL).hpp : generate_$(INTEGRAL).py
	python $<

$(INTEGRAL)/lib$(INTEGRAL).a : $(INTEGRAL)/$(INTEGRAL).hpp
	$(MAKE) -C $(INTEGRAL) static

$(INTEGRAL)/lib$(INTEGRAL)_pylink.so : $(INTEGRAL)/lib$(INTEGRAL).a
	$(MAKE) -C $(INTEGRAL) pylink
