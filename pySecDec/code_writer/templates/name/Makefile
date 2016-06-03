include Makefile.conf

compile : source $(LIBINTEGRANDS)

source:
	$(MAKE) -C codegen $@

test.o : $(LIBINTEGRANDS)

test : test.o $(LIBINTEGRANDS)
	$(CPP) $(CPPFLAGS) -o $@ $+ $(LDFLAGS)

very-clean : clean
	for dir in */; do $(MAKE) -C $$dir $@; done

clean:
	for dir in */; do $(MAKE) -C $$dir $@; done
	rm -f *.o test

$(LIBINTEGRANDS) : source
