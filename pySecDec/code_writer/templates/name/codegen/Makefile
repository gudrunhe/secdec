include ../Makefile.conf
INTEGRAND_STAMPS = $(patsubst sector%%.h,sector%%-stamp,$(wildcard sector*.h))
CONTOUR_DEFORMATION_STAMPS = $(patsubst contour_deformation_sector%%.h,contour_deformation_sector%%-stamp,$(wildcard contour_deformation_sector*.h))

compile : source

source : $(INTEGRAND_STAMPS) $(CONTOUR_DEFORMATION_STAMPS)

sector%%-stamp : sector%%.h
	@ # generate c++ code
	$(FORM) -D sectorID=$(patsubst sector%%.h,%%,$<) write_integrand.frm

	@ # move source code to appropriate directory
	mv sector_$(patsubst sector%%.h,%%,$<).hpp ../integrands
	mv sector_$(patsubst sector%%.h,%%,$<)_*.hpp ../integrands
	mv sector_$(patsubst sector%%.h,%%,$<)_*.cpp ../integrands

	@ # "mark" done
	touch $@

contour_deformation_sector%%-stamp : contour_deformation_sector%%.h
	@ # generate c++ code
	$(FORM) -D sectorID=$(patsubst contour_deformation_sector%%.h,%%,$<) write_contourdef.frm

	@ # move source code to appropriate directory
	mv contour_deformation_sector_$(patsubst contour_deformation_sector%%.h,%%,$<).hpp ../integrands
	mv contour_deformation_sector_$(patsubst contour_deformation_sector%%.h,%%,$<).cpp ../integrands

	@ # "mark" done
	touch $@

very-clean : clean
	rm -f *-stamp *.cpp *.hpp
