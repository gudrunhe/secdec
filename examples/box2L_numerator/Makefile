INTEGRALS = $(patsubst %.py,%,$(wildcard *.py))

# CXXFLAGS
ifeq (x$(SECDEC_CONTRIB), x)
CXXFLAGS =
else
CXXFLAGS = -I$(SECDEC_CONTRIB)/include
endif
CXXFLAGS += -c -std=c++11 -O2 -I.

# LDFLAGS
ifeq (x$(SECDEC_CONTRIB), x)
LDFLAGS =
else
LDFLAGS = -L$(SECDEC_CONTRIB)/lib
endif
LDFLAGS += $(patsubst %,-L%,$(INTEGRALS)) $(patsubst %,-l%,$(INTEGRALS)) -lgsl -lgslcblas -lcuba -lm

integrate_box2L : integrate_box2L.o $(patsubst %,lib%.a,$(INTEGRALS))
	$(CXX) -o $@ $< $(LDFLAGS)

integrate_box2L.o : integrate_box2L.cpp $(INTEGRALS)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(INTEGRALS) :
	python $@.py

lib%.a : %
	$(MAKE) -C $< $@

.PHONY : clean
clean :
	rm -f *.o *.so *.a integrate_box2L
	for dir in *; do if [ -d $$dir ]; then $(MAKE) -C $$dir $@; fi; done

.PHONY : very-clean
very-clean :
	rm -rf */ *.o *.so *.a integrate_box2L
