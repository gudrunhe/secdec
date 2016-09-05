include Makefile.conf

lib$(NAME).a lib$(NAME).so : source
	$(MAKE) -C src $@

source:
	$(MAKE) -C codegen $@

test.o : source

test : test.o lib$(NAME).a
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

very-clean : clean
	for dir in */; do $(MAKE) -C $$dir $@; done

clean:
	for dir in */; do $(MAKE) -C $$dir $@; done
	rm -f *.o *.so *.a test
