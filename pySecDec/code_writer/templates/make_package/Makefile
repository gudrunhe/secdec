ORIG_CXXFLAGS:=$(CXXFLAGS)

include Makefile.conf
include $(wildcard codegen/sector*.d)

source : $(SECTOR_CPP)

lib$(NAME).a : $(patsubst %%.cpp,%%.o,$(SECTOR_CPP)) src/integrands.o src/pole_structures.o src/prefactor.o
	ar rs $@ $+

lib$(NAME).so : $(patsubst %%.cpp,%%.o,$(SECTOR_CPP)) src/integrands.o src/pole_structures.o src/prefactor.o
	$(CXX) -o $@ -shared $+

QMC_TEMPLATE_OBJECTS = $(patsubst %%.cpp,%%.o,$(wildcard pylink/qmc_template_instantiations_*.cpp))

$(NAME)_pylink.so : pylink/pylink.o lib$(NAME).a $(QMC_TEMPLATE_OBJECTS)
	$(CXX) -shared -o $@ $+ $(LDFLAGS)

integrate_$(NAME) : integrate_$(NAME).o lib$(NAME).a
	$(CXX) -o $@ integrate_$(NAME).o lib$(NAME).a $(LDFLAGS)

ifdef SECDEC_WITH_CUDA_FLAGS
cuda_integrate_$(NAME) : cuda_integrate_$(NAME).o lib$(NAME).a
	$(CXX) -o $@ cuda_integrate_$(NAME).o lib$(NAME).a $(LDFLAGS)
endif

very-clean:: clean
	rm -f codegen/*.done src/*sector*.[ch]pp

clean::
	rm -f *.o *.so *.a pylink/*.o src/*.o integrate_$(NAME) cuda_integrate_$(NAME)
	rm -f distsrc/*.o distsrc/*.fatbin disteval/*.so disteval/*.fatbin

# implicit rule to build object files
%%.o : %%.cpp
ifdef SECDEC_WITH_CUDA_FLAGS
	$(CXX) -dc $(CXXFLAGS) -Xptxas -O0 -Xptxas --disable-optimizer-constants -Xcompiler -fPIC $< -o $@
else
	$(CXX) -c $(CXXFLAGS) -fPIC $< -o $@
endif

codegen/sector%%.done: codegen/sector%%.h
	@# generate c++ code
	cd codegen && $(PYTHON) formwrapper $(FORMCALL) -D sectorID=$(patsubst codegen/sector%%.h,%%,$<) write_integrand.frm
	$(PYTHON) codegen/export_sector $(patsubst %%.h,%%.info,$<) ./
	touch $@

# The following is for the distributed evaluation.

SECTOR_ORDERS:=$(patsubst src/sector_%%.cpp,%%,$(filter src/sector_%%.cpp,$(SECTOR_CPP)))
SECTOR_ORDERS:=$(foreach a,$(SECTOR_ORDERS),$(if $(findstring _,$a),$a,))

ifdef SECDEC_WITH_CUDA_FLAGS
disteval.done: disteval/$(NAME).fatbin disteval/builtin_cuda.fatbin disteval/$(NAME).so disteval/builtin_cpu.so
else
disteval.done: disteval/$(NAME).so disteval/builtin_cpu.so
endif
	date >$@

# CPU files (.so)

XCXXFLAGS=-std=c++14 -O3 -funsafe-math-optimizations $(ORIG_CXXFLAGS)

DIST_SO_OBJECTS = $(patsubst %%,distsrc/sector_%%.o,$(SECTOR_ORDERS))

distsrc/%%.o: distsrc/%%.cpp
	$(CXX) -c -o $@ -fPIC $(XCXXFLAGS) $^

disteval/$(NAME).so: $(DIST_SO_OBJECTS)
	$(CXX) -shared -o $@ $^

disteval/builtin_cpu.so: distsrc/builtin_cpu.o
	$(CXX) -shared -o $@ $^

# CUDA files (.fatbin)

NVCC?=nvcc
NVCCFLAGS=-std=c++14 $(SECDEC_WITH_CUDA_FLAGS)

DIST_FATBIN_OBJECTS = $(patsubst %%,distsrc/sector_%%.fatbin,$(SECTOR_ORDERS))

distsrc/%%.fatbin: distsrc/%%.cu
	$(NVCC) $(NVCCFLAGS) -dc -fatbin -o $@ $^

disteval/$(NAME).fatbin: $(DIST_FATBIN_OBJECTS)
	$(NVCC) $(NVCCFLAGS) -dlink -fatbin -o $@ $^

disteval/builtin_cuda.fatbin: distsrc/builtin_cuda.fatbin
	$(NVCC) $(NVCCFLAGS) -dlink -fatbin -o $@ $^
