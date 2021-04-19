INTEGRALS = $(patsubst %.py,%,$(wildcard *.py))

# CXXFLAGS
SECDEC_CONTRIB = $(shell python3 -m pySecDecContrib --dirname)
CXXFLAGS = -I"$(SECDEC_CONTRIB)/include" -c -std=c++11 -O2 -I.

# LDFLAGS
LDFLAGS = -L"$(SECDEC_CONTRIB)/lib" $(patsubst %,-L%,$(INTEGRALS)) $(patsubst %,-l%,$(INTEGRALS)) -lgsl -lgslcblas -lcuba -lm

amp : amp.o $(patsubst %,lib%.a,$(INTEGRALS))
	$(CXX) -o $@ $< $(LDFLAGS)

amp.o : amp.cpp $(INTEGRALS)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(INTEGRALS) :
	python3 $@.py

lib%.a : %
	$(MAKE) -C $< $@

.PHONY : clean
clean :
	rm -f *.o *.so *.a amp
	for dir in *; do if [ -d $$dir ]; then $(MAKE) -C $$dir $@; fi; done

.PHONY : very-clean
very-clean :
	rm -rf */ *.o *.so *.a amp
