include Makefile.conf

compile : source
	$(MAKE) -C integrands $@

source:
	$(MAKE) -C codegen $@

test.o : compile

test : test.o
	$(CPP) $(CPPFLAGS) -o $@ $+ $(TOPDIR)/integrands/integrands.a $(LDFLAGS)

very-clean : clean
	for dir in */; do $(MAKE) -C $$dir $@; done

clean:
	for dir in */; do $(MAKE) -C $$dir $@; done
	rm -f *.o test
