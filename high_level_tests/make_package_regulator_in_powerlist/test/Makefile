include ../Makefile.conf
include ../../test_suite/Makefile.conf

test.o : test.cpp $(TEST_FRAMEWORK_DIR)/catch_amalgamated.hpp
	$(CXX) $(CXXFLAGS) -I../$(INTEGRAL_NAME) $< -c -o $@

test.exe : test.o $(TEST_FRAMEWORK_DIR)/catch.o ../$(INTEGRAL_NAME)/lib$(INTEGRAL_NAME).a
	$(CXX) $+ -o $@ -L../$(INTEGRAL_NAME) -l$(INTEGRAL_NAME) $(LDFLAGS)

.PHONY : test
test : test.log
test.log : test.exe test.py
	(echo "CXX test:" && ./test.exe && echo "PYTHON test:" && $(PYTHON) test.py -v 2>&1) 2>&1 >test.log \
	&& printf "\n@@@ SUCCESS @@@" >> test.log || printf "\n@@@ FAILURE @@@" >> test.log

.PHONY : clean very-clean
clean :
	rm -rf *.o *.exe *.log

very-clean  : clean

