INTEGRALS = $(patsubst %.py,%,$(wildcard *.py))

# CXXFLAGS
SECDEC_CONTRIB = $(shell python3 -m pySecDecContrib --dirname)
CXXFLAGS = -I"$(SECDEC_CONTRIB)/include" -c -std=c++14 -O2 -I.

# LDFLAGS
LDFLAGS = -L"$(SECDEC_CONTRIB)/lib" $(patsubst %,-L%,$(INTEGRALS)) $(patsubst %,-l%,$(INTEGRALS)) -lgsl -lgslcblas -lcuba -lginac -lcln -lm

integrate_box2L : integrate_box2L.o $(patsubst %,lib%.a,$(INTEGRALS))
	$(CXX) -o $@ $< $(LDFLAGS)

integrate_box2L.o : integrate_box2L.cpp $(INTEGRALS)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(INTEGRALS) :
	python3 $@.py

lib%.a : %
	$(MAKE) -C $< $@

.PHONY : clean
clean :
	rm -f *.o *.so *.a integrate_box2L
	for dir in *; do if [ -d $$dir ]; then $(MAKE) -C $$dir $@; fi; done

.PHONY : very-clean
very-clean :
	rm -rf */ *.o *.so *.a integrate_box2L
