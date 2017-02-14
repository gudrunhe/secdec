include Makefile.conf

lib$(NAME).a lib$(NAME).so : source
	$(MAKE) -C src ../$@

source:
	$(MAKE) -C codegen $@

$(NAME).o : source

$(NAME) : $(NAME).o lib$(NAME).a
	$(CXX) $(CXXFLAGS) -o $@ $< -l$(NAME) $(LDFLAGS)

very-clean : clean
	for dir in */; do $(MAKE) -C $$dir $@; done

clean:
	for dir in */; do $(MAKE) -C $$dir $@; done
	rm -f *.o *.so *.a $(NAME)
