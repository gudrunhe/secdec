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

$(INTEGRALS_A): %%.a:
	$(MAKE) -C $(dir $@) $(basename $1)

$(NAME)_pylink.so: pylink/pylink.o src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A) $(QMC_TEMPLATE_OBJECTS)
	$(CXX) -shared -o $@ pylink/pylink.o src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A) $(QMC_TEMPLATE_OBJECTS) $(SUM_LDFLAGS)

lib$(NAME).a : src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A)
	( \
		echo "create $@"; \
		echo "addmod src/amplitude.o $(WINTEGRALS_OBJS)"; \
		for x in $(INTEGRALS_A); do echo "addlib $$x"; done; \
		echo "save"; \
		echo "end"; \
	) | ar -M

lib$(NAME).so : src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A)
ifdef SECDEC_WITH_CUDA_FLAGS
	$(CXX) -shared -o $@ src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A) $(SUM_LDFLAGS)
else
	$(CXX) -shared -o $@ src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A) $(SUM_LDFLAGS) -Wl,-undefined,dynamic_lookup
endif

# build the example executable
integrate_$(NAME) : integrate_$(NAME).o lib$(NAME).a $(INTEGRALS_A)
	$(CXX) -o $@ integrate_$(NAME).o lib$(NAME).a $(INTEGRALS_A) $(SUM_LDFLAGS)

very-clean :: clean
	for dir in */; do if [ -e "$$dir/Makefile" ]; then $(MAKE) -C "$$dir" $@; fi; done

clean ::
	for dir in */; do if [ -e "$$dir/Makefile" ]; then $(MAKE) -C "$$dir" $@; fi; done
	rm -f *.o *.so *.a pylink/*.o src/*.o integrate_$(NAME)
