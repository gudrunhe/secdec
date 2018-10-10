include Makefile.conf

$(foreach INTEGRAL,$(INTEGRALS),lib/lib$(INTEGRAL).so) lib/lib$(NAME).so:
	$(MAKE) -C src ../$@

source:
	for INTEGRAL in $(INTEGRALS); do $(MAKE) -C $${INTEGRAL} $@; done

# rule to build the integral libraries
libs:
	$(MAKE) -C src $@

# build the example executable
ifdef SECDEC_WITH_CUDA
integrate_$(NAME) : integrate_$(NAME).o $(foreach INTEGRAL,$(INTEGRALS),lib/lib$(INTEGRAL).so) lib/lib$(NAME).so
	$(CXX) -o $@ $< $(LDFLAGS) -Xlinker '-rpath lib' -l$(NAME) $(patsubst %%,-l%%,$(INTEGRALS))
else
integrate_$(NAME) : integrate_$(NAME).o $(foreach INTEGRAL,$(INTEGRALS),lib/lib$(INTEGRAL).so) lib/lib$(NAME).so
	$(CXX) -o $@ $< $(LDFLAGS) -Wl,-rpath,lib -l$(NAME) $(patsubst %%,-l%%,$(INTEGRALS))
endif

very-clean : clean
	for dir in */; do $(MAKE) -C $$dir $@; done

clean:
	for dir in */; do $(MAKE) -C $$dir $@; done
	rm -f *.o *.so *.a integrate_$(NAME)
