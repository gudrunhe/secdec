include ../Makefile.conf

OBJECTS = $(patsubst %%.cpp,%%.o,$(wildcard *.cpp))

../lib$(NAME).a : lib$(NAME).a
	mv $< $@

.INTERMEDIATE : lib$(NAME).a
lib$(NAME).a : $(OBJECTS)
	ar rs $@ $+

../lib$(NAME).so : lib$(NAME).so
	mv $< $@

.INTERMEDIATE : lib$(NAME).so
lib$(NAME).so : $(OBJECTS)
	$(CXX) -o $@ -shared $+ $(LDFLAGS)

clean:
	rm -f *.o lib$(NAME).a lib$(NAME).so

very-clean : clean
	rm -f *sector_*.cpp *sector_*.hpp

# implicit rule to build object files
ifdef SECDEC_WITH_CUDA
%%.o : %%.cpp
	$(CXX) -dc $(CXXFLAGS) -Xptxas -O0 -Xptxas --disable-optimizer-constants -Xcompiler -fPIC $< -o $@
else
%%.o : %%.cpp
	$(CXX) -c $(CXXFLAGS) -fPIC $< -o $@
endif
