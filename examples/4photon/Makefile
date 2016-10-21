include Makefile.conf

BUB=yyyy_bubble
BOX=yyyy_box6Dim

lib$(NAME_I2).a lib$(NAME_I2).so : $(BUB)/Makefile
	$(MAKE) -C $(BUB)/src $@
	ln $(BUB)/src/$@

lib$(NAME_I4).a lib$(NAME_I4).so : $(BOX)/Makefile
	$(MAKE) -C $(BOX)/src $@
	ln $(BOX)/src/$@

$(BUB)/Makefile :
	python 4y_bubble.py
	$(MAKE) -C $(BUB)/codegen source

$(BOX)/Makefile :
	python 4y_box.py
	$(MAKE) -C $(BOX)/codegen source

amp.o test.o : $(BUB)/Makefile $(BOX)/Makefile

test : amp.o lib$(NAME_I2).a lib$(NAME_I4).a
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

very-clean : clean
	for dir in */; do rm -rf $$dir $@; done

clean:
	for dir in */; do $(MAKE) -C $$dir $@; done
	rm -f *.o *.so *.a test
