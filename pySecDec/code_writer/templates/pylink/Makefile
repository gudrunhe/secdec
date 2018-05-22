include ../Makefile.conf

../$(NAME)_pylink.so : pylink.o ../lib$(NAME).a
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
