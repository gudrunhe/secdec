# get the absolute path to this directory
TOPDIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# common options for all packages
CONFIGURE_COMMON = CFLAGS=-fexceptions --prefix=$(TOPDIR)/install --libdir=$(TOPDIR)/lib --includedir=$(TOPDIR)/include --bindir=$(TOPDIR)/bin

# the included packages
CPP_PACKAGES = Cuba form gsl secdecutil


.PHONY : clean $(CPP_PACKAGES) pySecDec

.PHONY .SILENT : install
begin_green = '\033[0;32m'
end_green   = '\033[0m'
install : $(CPP_PACKAGES) pySecDec
	echo -e $(begin_green)"Successfully built \"pySecDec\" and its dependencies."$(end_green) ; \
	echo -e $(begin_green)"To use \"pySecDec\", add the following two commands to your \".bashrc\" or \".profile\":"$(end_green) ; \
	echo -e $(begin_green)"export SECDEC_CONTRIB=$(TOPDIR)"$(end_green) ; \
	echo -e $(begin_green)"export PYTHONPATH=$(TOPDIR)pySecDec-$(shell cd pySecDec-*; python -c 'import pySecDec; print(pySecDec.__version__)'):\$$PYTHONPATH"$(end_green)

pySecDec :
	tar -xf $@-*.tar.gz

CubaCONFIGURE = $(CONFIGURE_COMMON)
formCONFIGURE = $(CONFIGURE_COMMON)
gslCONFIGURE = $(CONFIGURE_COMMON) --disable-shared --enable-static
secdecutilCONFIGURE = $(CONFIGURE_COMMON)
$(CPP_PACKAGES) :
	tar -xf $@-*.tar.gz && \
	cd $@-* && \
	./configure $($@CONFIGURE) && \
	$(MAKE) && \
	$(MAKE) install

clean :
	rm -rf install lib include bin
	for PACKAGE_NAME in $(CPP_PACKAGES); do \
		rm -rf $${PACKAGE_NAME}-*/; \
	done

