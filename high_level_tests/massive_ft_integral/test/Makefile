include ../Makefile.conf
include ../../test_suite/Makefile.conf

ifdef SECDEC_WITH_CUDA
COMPILE = -dc
else
COMPILE = -c
endif

.PHONY : test
test : test.log
test.log : test.py
	(echo "PYTHON test:" && $(PYTHON) test.py -q) 2>&1 >test.log \
	&& printf "\n@@@ SUCCESS @@@" >> test.log || printf "\n@@@ FAILURE @@@" >> test.log

.PHONY : clean very-clean
clean :
	rm -rf *.o *.exe *.log

very-clean  : clean

