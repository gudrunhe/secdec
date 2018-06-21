include Makefile.conf

$(foreach INTEGRAL,$(INTEGRALS),lib/lib$(INTEGRAL).so) lib/lib$(NAME).so : config_$(NAME).hpp
	$(MAKE) -C src ../$@

source:
	for INTEGRAL in $(INTEGRALS); do $(MAKE) -C $${INTEGRAL} $@& done; wait

# rule to build the integral libraries
libs:
	$(MAKE) -C src $@

# build the example executable
ifdef SECDEC_WITH_CUDA
ifneq (x$(SECDEC_CONTRIB), x)
CONTRIB_RPATH := -Xlinker "-rpath $(SECDEC_CONTRIB)/lib"
else
CONTRIB_RPATH :=
endif
integrate_$(NAME) : integrate_$(NAME).o $(foreach INTEGRAL,$(INTEGRALS),lib/lib$(INTEGRAL).so) lib/lib$(NAME).so
	$(CXX) -o $@ $< $(LDFLAGS) -Xlinker '-rpath lib' $(CONTRIB_RPATH) -l$(NAME) $(patsubst %%,-l%%,$(INTEGRALS))
else
ifneq (x$(SECDEC_CONTRIB), x)
CONTRIB_RPATH := -Wl,-rpath,$(SECDEC_CONTRIB)/lib
else
CONTRIB_RPATH :=
endif
integrate_$(NAME) : integrate_$(NAME).o $(foreach INTEGRAL,$(INTEGRALS),lib/lib$(INTEGRAL).so) lib/lib$(NAME).so
	$(CXX) -o $@ $< $(LDFLAGS) -Wl,-rpath,lib $(CONTRIB_RPATH) -l$(NAME) $(patsubst %%,-l%%,$(INTEGRALS))
endif

very-clean : clean
	for dir in */; do $(MAKE) -C $$dir $@; done

clean:
	for dir in */; do $(MAKE) -C $$dir $@; done
	rm -f *.o *.so *.a integrate_$(NAME)