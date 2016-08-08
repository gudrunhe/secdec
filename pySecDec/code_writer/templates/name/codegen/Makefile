include ../Makefile.conf
INTEGRAND_STAMPS = $(patsubst sector%%.h,sector%%-stamp,$(wildcard sector*.h))

compile : source

source : $(INTEGRAND_STAMPS)

sector%%-stamp : sector%%.h
	@ # generate c++ code
	$(FORM) -D sectorID=$(patsubst sector%%.h,%%,$<) write_integrand.frm

	@ # move source code to appropriate directory
	$(PYTHON) mv.py sector_$(patsubst sector%%.h,%%,$<).hpp ../integrands
	$(PYTHON) mv.py *sector_$(patsubst sector%%.h,%%,$<)_*.hpp ../integrands
	$(PYTHON) mv.py *sector_$(patsubst sector%%.h,%%,$<)_*.cpp ../integrands

	@ # "mark" done
	touch $@

very-clean : clean
	rm -f *-stamp *.cpp *.hpp
