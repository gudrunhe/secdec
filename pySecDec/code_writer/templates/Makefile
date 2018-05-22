include Makefile.conf

lib$(NAME).a lib$(NAME).so : source
	$(MAKE) -C src ../$@

$(NAME)_pylink.so : lib$(NAME).a
	$(MAKE) -C pylink ../$@

source:
	$(MAKE) -C codegen $@

integrate_$(NAME).o : source

integrate_$(NAME) : integrate_$(NAME).o lib$(NAME).a
	$(CXX) -O2 $< -o $@ -l$(NAME) $(LDFLAGS)

ifdef SECDEC_WITH_CUDA
cuda_integrate_$(NAME) : cuda_integrate_$(NAME).o lib$(NAME).a
	$(CXX) -O2 $< -o $@ -l$(NAME) $(LDFLAGS)
endif

very-clean : clean
	for dir in */; do $(MAKE) -C $$dir $@; done

clean:
	for dir in */; do $(MAKE) -C $$dir $@; done
	rm -f *.o *.so *.a integrate_$(NAME) cuda_integrate_$(NAME)

# implicit rule to build object files
ifdef SECDEC_WITH_CUDA
%%.o : %%.cpp
	$(CXX) -dc $(CXXFLAGS) -Xcompiler -fPIC $< -o $@
else
%%.o : %%.cpp
	$(CXX) -c $(CXXFLAGS) -fPIC $< -o $@
endif
