# get the absolute path to this directory
TOPDIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# common options for all packages
CONFIGURE_COMMON = CFLAGS=-fexceptions --prefix=$(TOPDIR)/install --libdir=$(TOPDIR)/lib --includedir=$(TOPDIR)/include --bindir=$(TOPDIR)/bin

# the included packages
PACKAGES = Cuba form gsl secdecutil


.PHONY : all clean $(PACKAGES)

all : $(PACKAGES)


CubaVERSION = 4.2
CubaCONFIGURE = $(CONFIGURE_COMMON)

formVERSION = 4.1
formCONFIGURE = $(CONFIGURE_COMMON)

gslVERSION = 2.1
gslCONFIGURE = $(CONFIGURE_COMMON) --disable-shared --enable-static

secdecutilVERSION = 0.1.1
secdecutilCONFIGURE = $(CONFIGURE_COMMON)


$(PACKAGES) :
	tar -xf $@-$($@VERSION).tar.gz && \
	cd $@-$($@VERSION) && \
	./configure $($@CONFIGURE) && \
	cd .. && \
	$(MAKE) -C $@-$($@VERSION) && \
	$(MAKE) -C $@-$($@VERSION) install

# export all variables to the shell
export
clean :
	rm -rf install lib include bin
	for PACKAGE_NAME in $(PACKAGES); do \
		eval 'VERSION_NUMBER=$${PACKAGE_NAME}VERSION' && \
		eval 'VERSION_NUMBER=$${'$$VERSION_NUMBER'}' && \
		rm -rf $${PACKAGE_NAME}-$${VERSION_NUMBER}; \
	done

