include Makefile.conf

WINTEGRALS_OBJS = $(foreach INTEGRAL,$(INTEGRALS),src/$(INTEGRAL)_weighted_integral.o)
INTEGRALS_SO = $(foreach INTEGRAL,$(INTEGRALS),lib/lib$(INTEGRAL).so)
INTEGRALS_A = $(foreach INTEGRAL,$(INTEGRALS),$(INTEGRAL)/lib$(INTEGRAL).a)

source:
	for INTEGRAL in $(INTEGRALS); do $(MAKE) -C $${INTEGRAL} $@& done; wait

# rule to build the integral libraries
libs:
	$(MAKE) -C src $@

pylink: libs $(NAME)_pylink.so


$(INTEGRALS_A) :
	for INTEGRAL in $(INTEGRALS); do $(MAKE) -C $${INTEGRAL} static& done; wait

# build the python shared library
$(NAME)_pylink.so: pylink/pylink.o src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A)
	$(CXX) -shared -o $@ pylink/pylink.o src/amplitude.o $(WINTEGRALS_OBJS) $(foreach INTEGRAL,$(INTEGRALS),$(INTEGRAL)/src/*.o) $(LDFLAGS)

# build the example executable
ifdef SECDEC_WITH_CUDA_FLAGS
ifneq (x$(SECDEC_CONTRIB), x)
CONTRIB_RPATH := -Xlinker "-rpath $(SECDEC_CONTRIB)/lib"
else
CONTRIB_RPATH :=
endif
integrate_$(NAME) : integrate_$(NAME).o $(INTEGRALS_SO) lib/lib$(NAME).so
	$(CXX) -o $@ $< $(LDFLAGS) -Xlinker '-rpath lib' $(CONTRIB_RPATH) -l$(NAME) $(patsubst %%,-l%%,$(INTEGRALS))
else
ifneq (x$(SECDEC_CONTRIB), x)
CONTRIB_RPATH := -Wl,-rpath,$(SECDEC_CONTRIB)/lib
else
CONTRIB_RPATH :=
endif
integrate_$(NAME) : integrate_$(NAME).o $(INTEGRALS_SO) lib/lib$(NAME).so
	$(CXX) -o $@ $< $(LDFLAGS) -Wl,-rpath,lib $(CONTRIB_RPATH) -l$(NAME) $(patsubst %%,-l%%,$(INTEGRALS))
endif

very-clean : clean
	for dir in */; do $(MAKE) -C $$dir $@; done

clean:
	for dir in */; do $(MAKE) -C $$dir $@; done
	rm -f *.o *.so *.a integrate_$(NAME)
