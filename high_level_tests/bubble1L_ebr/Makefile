include ./Makefile.conf
include ../test_suite/Makefile.conf

$(INTEGRAL_NAME) : $(INTEGRAL_NAME).py
	@rm -rf $(INTEGRAL_NAME)
	$(PYTHON) $<

.PHONY : static-library
static-library : $(INTEGRAL_NAME)
	$(MAKE) -C $(INTEGRAL_NAME) static

.PHONY : pylink-library
pylink-library : $(INTEGRAL_NAME) static-library
	$(MAKE) -C $(INTEGRAL_NAME) pylink

.PHONY : test
test : static-library pylink-library
	echo $(EXAMPLES_DIR)
	$(MAKE) -C test INTEGRAL_NAME=$(INTEGRAL_NAME) $@

.PHONY : clean
clean :
	$(MAKE) -C $(INTEGRAL_NAME) $@
	$(MAKE) -C test INTEGRAL_NAME=$(INTEGRAL_NAME) $@

.PHONY : very-clean
very-clean :
	rm -rf $(INTEGRAL_NAME)
	$(MAKE) -C test INTEGRAL_NAME=$(INTEGRAL_NAME) $@

