include ../Makefile.conf

OBJECTS = $(patsubst %%.cpp,%%.o,$(wildcard *.cpp))

compile : source integrands.a

# TODO: add option for dynamic library
integrands.a : $(OBJECTS)
	ar rs $@ $+

clean:
	rm -f *.o

very-clean : clean
	rm -f sector_*.cpp sector_*.hpp
