# get the absolute path to this directory
TOPDIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# common options for all packages
CONFIGURE_COMMON = --prefix=$(TOPDIR)/install --libdir=$(TOPDIR)/lib --includedir=$(TOPDIR)/include --bindir=$(TOPDIR)/bin

# the included packages
CPP_PACKAGES = Cuba form gsl secdecutil


.PHONY : clean $(CPP_PACKAGES) pySecDec nauty

.PHONY .SILENT : install
begin_green = '\033[0;32m'
end_green   = '\033[0m'
install : nauty $(CPP_PACKAGES) pySecDec
	printf $(begin_green)"Successfully built \"pySecDec\" and its dependencies.\n"$(end_green) ; \
	printf $(begin_green)"To use \"pySecDec\", add the following two commands to your \".bashrc\" or \".profile\":\n"$(end_green) ; \
	printf $(begin_green)"export SECDEC_CONTRIB=$(TOPDIR)\n"$(end_green) ; \
	printf $(begin_green)"export PYTHONPATH=$(TOPDIR)pySecDec-$(shell cd pySecDec-*/; python -c 'import pySecDec; print(pySecDec.__version__)'):\$$PYTHONPATH\n"$(end_green)

pySecDec :
	tar -xf $@-*.tar.gz

CubaCONFIGURE = $(CONFIGURE_COMMON) CFLAGS="-O3 -fomit-frame-pointer -ffast-math -fno-finite-math-only -fexceptions -fPIC"
formCONFIGURE = $(CONFIGURE_COMMON)
gslCONFIGURE = $(CONFIGURE_COMMON) --disable-shared --enable-static --with-pic
secdecutilCONFIGURE = $(CONFIGURE_COMMON)
$(CPP_PACKAGES) :
	tar -xf $@-*.tar.gz && \
	cd $@-*/ && \
	./configure $($@CONFIGURE) && \
	$(MAKE) && \
	$(MAKE) install

nauty :
	tar -xf $@*.tar.gz && \
	cd $@*/ && \
	./configure && \
	$(MAKE) dreadnaut && \
	mkdir -p $(TOPDIR)/bin/ && \
	cp dreadnaut $(TOPDIR)/bin/
clean :
	rm -rf install lib include bin nauty*/
	for PACKAGE_NAME in $(CPP_PACKAGES); do \
		rm -rf $${PACKAGE_NAME}-*/; \
	done

