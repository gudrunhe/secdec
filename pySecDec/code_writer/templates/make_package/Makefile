include Makefile.conf
include $(wildcard codegen/sector*.d)

source : $(SECTOR_CPP)

lib$(NAME).a : $(patsubst %%.cpp,%%.o,$(SECTOR_CPP)) src/integrands.o src/pole_structures.o src/prefactor.o
	@rm -f $@
	lib=$$(mktemp) && \
		rm -f "$$lib" && \
		echo src/*.o | xargs $(AR) -c -q "$$lib" && \
		$(AR) -s "$$lib" && \
		mv "$$lib" $@

lib$(NAME).so : lib$(NAME).a
	$(XCC) -o $@ -shared $+ $(XLDFLAGS)

QMC_TEMPLATE_OBJECTS = $(patsubst %%.cpp,%%.o,$(wildcard pylink/qmc_template_instantiations_*.cpp))

$(NAME)_pylink.so : pylink/pylink.o lib$(NAME).a $(QMC_TEMPLATE_OBJECTS)
	$(XCC) -shared -o $@ $+ $(XLDFLAGS)

integrate_$(NAME) : integrate_$(NAME).o lib$(NAME).a
	$(XCC) -o $@ integrate_$(NAME).o lib$(NAME).a $(XLDFLAGS)

ifdef SECDEC_WITH_CUDA_FLAGS
cuda_integrate_$(NAME) : cuda_integrate_$(NAME).o lib$(NAME).a
	$(XCC) -o $@ cuda_integrate_$(NAME).o lib$(NAME).a $(XLDFLAGS)
endif

very-clean:: clean
	rm -f codegen/*.done src/*sector*.[ch]pp

clean::
	rm -f *.o *.so *.a pylink/*.o src/*.o integrate_$(NAME) cuda_integrate_$(NAME)
	rm -f disteval.done distsrc/*.o distsrc/*.fatbin disteval/*.so disteval/*.fatbin

# implicit rule to build object files
%%.o : %%.cpp
ifdef SECDEC_WITH_CUDA_FLAGS
	$(XCC) -dc $(XCCFLAGS) -Xptxas "-O0 --disable-optimizer-constants" -Xcompiler -fPIC $< -o $@
else
	$(XCC) -c $(XCCFLAGS) -fPIC $< -o $@
endif

codegen/sector%%.done: codegen/sector%%.h
	@# generate c++ code
	cd codegen && $(PYTHON) '$(SECDEC_CONTRIB)/bin/formwrapper' $(FORMCALL) -D sectorID=$(patsubst codegen/sector%%.h,%%,$<) '$(SECDEC_CONTRIB)/lib/write_integrand.frm'
	$(PYTHON) '$(SECDEC_CONTRIB)/bin/export_sector' $(patsubst %%.h,%%.info,$<) ./
	touch $@

# The following is for the distributed evaluation.

SECTOR_ORDERS:=$(patsubst src/sector_%%.cpp,%%,$(filter src/sector_%%.cpp,$(SECTOR_CPP)))
SECTOR_ORDERS:=$(foreach a,$(SECTOR_ORDERS),$(if $(findstring _,$a),$a,))

disteval: disteval.done

ifdef SECDEC_WITH_CUDA_FLAGS
disteval.done: disteval/$(NAME).fatbin disteval/builtin.fatbin disteval/$(NAME).so disteval/builtin.so
else
disteval.done: disteval/$(NAME).so disteval/builtin.so
endif
	date >$@

# CPU files (.so)

XCXXFLAGS=-std=c++14 -O3 -funsafe-math-optimizations $(CXXFLAGS)

DIST_SO_OBJECTS = $(patsubst %%,distsrc/sector_%%.o,$(SECTOR_ORDERS))

distsrc/%%.o: distsrc/%%.cpp
	$(CXX) -c -o $@ -fPIC $(XCXXFLAGS) $^

disteval/$(NAME).so: $(DIST_SO_OBJECTS)
	@echo distsrc/sector_*.o >$@.sourcelist
	$(CXX) -shared -o $@ @$@.sourcelist
	@rm -f $@.sourcelist

disteval/builtin.so: distsrc/builtin.o
	$(CXX) -shared -o $@ $^

# CUDA files (.fatbin)

XNVCCFLAGS=-std=c++14 $(SECDEC_WITH_CUDA_FLAGS) $(NVCCFLAGS)

DIST_FATBIN_OBJECTS = $(patsubst %%,distsrc/sector_%%.fatbin,$(SECTOR_ORDERS))

distsrc/%%.fatbin: distsrc/%%.cu
	$(NVCC) $(XNVCCFLAGS) -dc -fatbin -o $@ $^

disteval/$(NAME).fatbin: $(DIST_FATBIN_OBJECTS)
	@echo distsrc/sector_*.fatbin >$@.sourcelist
	$(NVCC) $(XNVCCFLAGS) -dlink -fatbin -o $@ --options-file $@.sourcelist
	@rm -f $@.sourcelist

disteval/builtin.fatbin: distsrc/builtin.fatbin
	$(NVCC) $(XNVCCFLAGS) -dlink -fatbin -o $@ $^
