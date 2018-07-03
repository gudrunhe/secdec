include ../Makefile.conf

../$(NAME)_pylink.so : pylink.o ../lib$(NAME).a qmc_template_instantiations1.o \
                       qmc_template_instantiations2.o qmc_template_instantiations3.o qmc_template_instantiations4.o \
                       qmc_template_instantiations5.o qmc_template_instantiations6.o qmc_template_instantiations7.o
	$(CXX) -shared -o $@ $+ $(LDFLAGS)

clean:
	rm -f *.o lib$(NAME).a lib$(NAME).so

# implicit rule to build object files
ifdef SECDEC_WITH_CUDA
%%.o : %%.cpp
	$(CXX) -dc $(CXXFLAGS) -Xcompiler -fPIC $< -o $@
else
%%.o : %%.cpp
	$(CXX) -c $(CXXFLAGS) -fPIC $< -o $@
endif
