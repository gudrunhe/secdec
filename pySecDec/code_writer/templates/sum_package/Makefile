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
	$(XCC) -shared -o $@ pylink/pylink.o src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A) $(QMC_TEMPLATE_OBJECTS) $(XLDFLAGS)

lib$(NAME).a : src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A)
	@rm -f $@
	dir=$$(mktemp -d) && \
		$(AR) -c -q "$$dir/lib.ar" src/amplitude.o $(WINTEGRALS_OBJS) && \
		cd "$$dir" && \
		$(foreach A,$(INTEGRALS_A),\
			$(AR) -x "$(CURDIR)/$(A)" && \
			echo *.o | xargs $(AR) -q lib.ar && \
			rm -f *.o && ) \
		$(AR) -s lib.ar && \
		mv lib.ar "$(CURDIR)/$@" && \
		rm -rf "$$dir"

lib$(NAME).so : src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A)
ifdef SECDEC_WITH_CUDA_FLAGS
	$(XCC) -shared -o $@ src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A) $(XLDFLAGS)
else
	$(XCC) -shared -o $@ src/amplitude.o $(WINTEGRALS_OBJS) $(INTEGRALS_A) $(XLDFLAGS) -Wl,-undefined,dynamic_lookup
endif

# build the example executable
integrate_$(NAME) : integrate_$(NAME).o lib$(NAME).a
	$(XCC) -o $@ integrate_$(NAME).o lib$(NAME).a $(XLDFLAGS)

very-clean :: clean
	for dir in */; do if [ -e "$$dir/Makefile" ]; then $(MAKE) -C "$$dir" $@; fi; done

clean ::
	for dir in */; do if [ -e "$$dir/Makefile" ]; then $(MAKE) -C "$$dir" $@; fi; done
	rm -f *.o *.so *.a pylink/*.o src/*.o integrate_$(NAME)
	rm -f disteval.done disteval/*.so disteval/*.fatbin $(foreach I,$(INTEGRALS),disteval/$I.json)

# implicit rule to build object files
ifdef SECDEC_WITH_CUDA_FLAGS
%%.o : %%.cpp
	$(XCC) -dc $(XCCFLAGS) -Xcompiler -fPIC $< -o $@
else
%%.o : %%.cpp
	$(XCC) -c $(XCCFLAGS) -fPIC $< -o $@
endif

# The following is for the distributed evaluation.

$(foreach I,$(INTEGRALS),disteval/$I.json): disteval/%%.json:; ln -f $*/$@ $@
$(foreach I,$(INTEGRALS),disteval/$I.so): disteval/%%.so: %%/disteval.done; ln -f $*/$@ $@
$(foreach I,$(INTEGRALS),disteval/$I.fatbin): disteval/%%.fatbin: %%/disteval.done; ln -f $*/$@ $@

disteval/builtin.so: $(word 1,$(INTEGRALS))/disteval.done; ln -f $(word 1,$(INTEGRALS))/$@ $@
disteval/builtin.fatbin: $(word 1,$(INTEGRALS))/disteval.done; ln -f $(word 1,$(INTEGRALS))/$@ $@

$(foreach I,$(INTEGRALS),$I/disteval.done)::
	$(MAKE) -C $(dir $@) disteval.done

disteval: disteval.done

ifdef SECDEC_WITH_CUDA_FLAGS
disteval.done: $(foreach I,$(INTEGRALS),disteval/$I.json disteval/$I.so disteval/$I.fatbin) disteval/builtin.so disteval/builtin.fatbin
else
disteval.done: $(foreach I,$(INTEGRALS),disteval/$I.json disteval/$I.so) disteval/builtin.so
endif
	date >$@
