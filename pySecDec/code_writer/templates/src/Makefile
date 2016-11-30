include ../Makefile.conf

OBJECTS = $(patsubst %%.cpp,%%.o,$(wildcard *.cpp))

../lib$(NAME).a : $(OBJECTS)
	ar rs $@ $+

../lib$(NAME).so : $(OBJECTS)
	$(CXX) -o $@ -shared $+

clean:
	rm -f *.o lib$(NAME).a lib$(NAME).so

very-clean : clean
	rm -f *sector_*.cpp *sector_*.hpp
