# get the absolute path to this directory
TOPDIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# common options for all packages
CONFIGURE_COMMON = --prefix=$(TOPDIR)/install --libdir=$(TOPDIR)/lib --includedir=$(TOPDIR)/include --bindir=$(TOPDIR)/bin

# the included autotools packages
AUTOTOOLS_PACKAGES = Cuba form gsl secdecutil

# the distributed version of pySecDec
PYSECDEC_VERSION = $(patsubst pySecDec-%.tar.gz,%,$(wildcard pySecDec-*.tar.gz))

.PHONY : clean $(AUTOTOOLS_PACKAGES) pySecDec nauty

.PHONY .SILENT : install
begin_green = '\033[0;32m'
end_green   = '\033[0m'
install : nauty $(AUTOTOOLS_PACKAGES) pySecDec
	echo -e $(begin_green)"Successfully built \"pySecDec\" and its dependencies."$(end_green) ; \
	echo -e $(begin_green)"To use \"pySecDec\", add the following two commands to your \".bashrc\" or \".profile\":"$(end_green) ; \
	echo -e $(begin_green)"export SECDEC_CONTRIB=$(TOPDIR)"$(end_green) ; \
	echo -e $(begin_green)"export PYTHONPATH=$(TOPDIR)pySecDec-$(PYSECDEC_VERSION):\$$PYTHONPATH"$(end_green)

pySecDec :
	tar -xf $@-*.tar.gz

CubaCONFIGURE = $(CONFIGURE_COMMON) CFLAGS="-O3 -fomit-frame-pointer -ffast-math -fno-finite-math-only -fexceptions -fPIC"
formCONFIGURE = $(CONFIGURE_COMMON)
gslCONFIGURE = $(CONFIGURE_COMMON) --disable-shared --enable-static --with-pic
secdecutilCONFIGURE = $(CONFIGURE_COMMON)

$(AUTOTOOLS_PACKAGES) :
	tar -xf $@-*.tar.gz && \
	cd $@-*/ && \
	./configure $($@CONFIGURE) && \
	(if [ "$@" == "Cuba" ]; then $(MAKE) -j1; else $(MAKE); fi) && \
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
	for PACKAGE_NAME in $(AUTOTOOLS_PACKAGES); do \
		rm -rf $${PACKAGE_NAME}-*/; \
	done

