include Makefile.conf

WINTEGRALS_OBJS = $(foreach INTEGRAL,$(INTEGRALS),src/$(INTEGRAL)_weighted_integral.o)
INTEGRALS_SO = $(foreach INTEGRAL,$(INTEGRALS),lib/lib$(INTEGRAL).so)
INTEGRALS_A = $(foreach INTEGRAL,$(INTEGRALS),$(INTEGRAL)/lib$(INTEGRAL).a)
QMC_TEMPLATE_OBJECTS = $(patsubst %%.cpp,%%.o,$(wildcard pylink/qmc_template_instantiations_*.cpp))

source:
	for INTEGRAL in $(INTEGRALS); do $(MAKE) -C $${INTEGRAL} $@& done; wait

# rule to build the integral libraries
libs:
	$(MAKE) -C src $@

pylink: libs $(NAME)_pylink.so


$(INTEGRALS_A) :
	for INTEGRAL in $(INTEGRALS); do $(MAKE) -C $${INTEGRAL} static& done; wait

# build the python shared library
$(NAME)_pylink.so: pylink/pylink.o src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A) $(QMC_TEMPLATE_OBJECTS)
	$(CXX) -shared -o $@ pylink/pylink.o src/amplitude.o $(WINTEGRALS_OBJS) $(foreach INTEGRAL,$(INTEGRALS),$(INTEGRAL)/src/*.o) $(LDFLAGS)

# build the example executable
ifdef SECDEC_WITH_CUDA_FLAGS
CONTRIB_RPATH := -Xlinker "-rpath $(SECDEC_CONTRIB)/lib"
integrate_$(NAME) : integrate_$(NAME).o $(INTEGRALS_SO) lib/lib$(NAME).so
	$(CXX) -o $@ $< $(LDFLAGS) -Xlinker '-rpath lib' $(CONTRIB_RPATH) -l$(NAME) $(patsubst %%,-l%%,$(INTEGRALS))
else
CONTRIB_RPATH := -Wl,-rpath,$(SECDEC_CONTRIB)/lib
integrate_$(NAME) : integrate_$(NAME).o $(INTEGRALS_SO) lib/lib$(NAME).so
	$(CXX) -o $@ $< $(LDFLAGS) -Wl,-rpath,lib $(CONTRIB_RPATH) -l$(NAME) $(patsubst %%,-l%%,$(INTEGRALS))
endif

very-clean : clean
	for dir in */; do $(MAKE) -C $$dir $@; done

clean:
	for dir in */; do $(MAKE) -C $$dir $@; done
	rm -f *.o *.so *.a integrate_$(NAME)
