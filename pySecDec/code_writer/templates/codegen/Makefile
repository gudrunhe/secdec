include ../Makefile.conf
INTEGRAND_STAMPS = $(patsubst sector%%.h,sector%%-stamp,$(wildcard sector*.h))

source : $(INTEGRAND_STAMPS)

sector%%-stamp : sector%%.h
	@ # generate c++ code
	$(FORMCALL) -D sectorID=$(patsubst sector%%.h,%%,$<) write_integrand.frm

	@ # move source code to appropriate directory
	$(PYTHON) mv.py sector_$(patsubst sector%%.h,%%,$<).cpp ../src
	$(PYTHON) mv.py *sector_$(patsubst sector%%.h,%%,$<)_*.hpp ../src
	$(PYTHON) mv.py *sector_$(patsubst sector%%.h,%%,$<)_*.cpp ../src

	@ # "mark" done
	touch $@

very-clean : clean
	rm -f *-stamp *.cpp *.hpp
