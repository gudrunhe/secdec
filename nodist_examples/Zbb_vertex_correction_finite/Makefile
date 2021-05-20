INTEGRALS = $(patsubst %.py,%,$(wildcard *.py))

# CXXFLAGS
SECDEC_CONTRIB = $(shell python3 -m pySecDecContrib --dirname)
CXXFLAGS = -I"$(SECDEC_CONTRIB)/include" -c -std=c++14 -O2 -I.

# LDFLAGS
LDFLAGS = -L"$(SECDEC_CONTRIB)/lib" $(patsubst %,-L%,$(INTEGRALS)) $(patsubst %,-l%,$(INTEGRALS)) -lcuba -lm

finite : finite.o $(patsubst %,lib%.a,$(INTEGRALS))
	$(CXX) -o $@ $< $(LDFLAGS)

finite.o : finite.cpp $(INTEGRALS)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(INTEGRALS) :
	python3 $@.py

lib%.a : %
	$(MAKE) -C $< $@

very-clean : clean
	for dir in */; do rm -rf $$dir $@; done

clean:
	for dir in */; do $(MAKE) -C $$dir $@; done
	rm -f *.o *.so *.a finite
