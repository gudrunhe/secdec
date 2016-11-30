include Makefile.conf

lib$(NAME).a lib$(NAME).so : source
	$(MAKE) -C src ../$@

source:
	$(MAKE) -C codegen $@

integrate_$(NAME).o : source

integrate_$(NAME) : integrate_$(NAME).o lib$(NAME).a
	$(CXX) $(CXXFLAGS) -o $@ $< -l$(NAME) $(LDFLAGS)

very-clean : clean
	for dir in */; do $(MAKE) -C $$dir $@; done

clean:
	for dir in */; do $(MAKE) -C $$dir $@; done
	rm -f *.o *.so *.a integrate_$(NAME)
