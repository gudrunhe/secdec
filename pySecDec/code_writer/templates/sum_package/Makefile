include Makefile.conf

WINTEGRALS_OBJS = $(foreach INTEGRAL,$(INTEGRALS),src/$(INTEGRAL)_weighted_integral.o)
INTEGRALS_A = $(foreach INTEGRAL,$(INTEGRALS),$(INTEGRAL)/lib$(INTEGRAL).a)
QMC_TEMPLATE_OBJECTS = $(patsubst %%.cpp,%%.o,$(wildcard pylink/qmc_template_instantiations_*.cpp))

# alias for the python shared library
pylink: $(NAME)_pylink.so

# alias for the static library
static: lib$(NAME).a

# alias for the dynamic library
dynamic: lib$(NAME).so

$(INTEGRALS_A):
	$(MAKE) -C $(dir $@) $(notdir $@)

$(NAME)_pylink.so: pylink/pylink.o src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A) $(QMC_TEMPLATE_OBJECTS)
	$(CXX) -shared -o $@ pylink/pylink.o src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A) $(QMC_TEMPLATE_OBJECTS) $(SUM_LDFLAGS)

lib$(NAME).a : src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A)
	@rm -f $@
	ar -c -q $@ src/amplitude.o $(WINTEGRALS_OBJS)
	dir=$$(mktemp -d) && \
	cd $$dir && \
	$(foreach A,$(INTEGRALS_A),ar -x "$(CURDIR)/$(A)" && ar -q "$(CURDIR)/$@" *.o && rm -f *.o && ) \
	rm -rf $$dir

lib$(NAME).so : src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A)
ifdef SECDEC_WITH_CUDA_FLAGS
	$(CXX) -shared -o $@ src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A) $(SUM_LDFLAGS)
else
	$(CXX) -shared -o $@ src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A) $(SUM_LDFLAGS) -Wl,-undefined,dynamic_lookup
endif

# build the example executable
integrate_$(NAME) : integrate_$(NAME).o lib$(NAME).a
	$(CXX) -o $@ integrate_$(NAME).o lib$(NAME).a $(SUM_LDFLAGS)

very-clean :: clean
	for dir in */; do if [ -e "$$dir/Makefile" ]; then $(MAKE) -C "$$dir" $@; fi; done

clean ::
	for dir in */; do if [ -e "$$dir/Makefile" ]; then $(MAKE) -C "$$dir" $@; fi; done
	rm -f *.o *.so *.a pylink/*.o src/*.o integrate_$(NAME)
