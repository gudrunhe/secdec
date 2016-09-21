include ../Makefile.conf

STATIC_OBJECTS = $(patsubst %%.cpp,%%.static.o,$(wildcard *.cpp))
DYNAMIC_OBJECTS = $(patsubst %%.cpp,%%.dynamic.o,$(wildcard *.cpp))

lib$(NAME).a : $(STATIC_OBJECTS)
	ar rs $@ $+
	cd .. && ln src/$@

lib$(NAME).so : $(DYNAMIC_OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ -shared $+
	cd .. && ln src/$@

clean:
	rm -f *.o lib$(NAME).a lib$(NAME).so

very-clean : clean
	rm -f *sector_*.cpp *sector_*.hpp
