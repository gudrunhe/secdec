include ./Makefile.conf

.DEFAULT_GOAL = catch.o

catch.o : catch.cpp catch.hpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

.PHONY : clean very-clean
clean :
	rm -rf *.o

very-clean  : clean
