include ../Makefile.conf

libs : ../lib/lib$(NAME).so $(foreach INTEGRAL,$(INTEGRALS),../lib/lib$(INTEGRAL).so)

.INTERMEDIATE : $(foreach INTEGRAL,$(INTEGRALS),lib$(INTEGRAL).so) lib$(NAME).so
$(foreach INTEGRAL,$(INTEGRALS),lib/lib$(INTEGRAL).so) : lib/lib%%.so : %%_weighted_integral.o
	export INTEGRAL=$(patsubst lib/lib%%.so,%%,$@) && \
	$(MAKE) -C ../$${INTEGRAL} static && \
	$(CXX) -o $@ -shared $${INTEGRAL}_weighted_integral.o ../$${INTEGRAL}/src/*.o $(LDFLAGS)
lib/lib$(NAME).so : amplitude.o
	$(CXX) -o $@ -shared $+ $(LDFLAGS) -Wl,-undefined,dynamic_lookup

$(foreach INTEGRAL,$(INTEGRALS),../lib/lib$(INTEGRAL).so) ../lib/lib$(NAME).so : ../lib/%% : lib/%%
	mv $< $@

clean :
	rm -f *.o *.so *.a

very-clean : clean
