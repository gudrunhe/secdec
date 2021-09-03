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
	$(PYTHON) codegen/export_sector $(patsubst %%.h,%%.info,$<) src/
	touch $@
