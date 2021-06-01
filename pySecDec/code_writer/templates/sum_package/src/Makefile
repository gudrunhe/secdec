include ../Makefile.conf

libs : ../lib/lib$(NAME).so $(foreach INTEGRAL,$(INTEGRALS),../lib/lib$(INTEGRAL).so)

.INTERMEDIATE : $(foreach INTEGRAL,$(INTEGRALS),lib$(INTEGRAL).so) lib$(NAME).so
$(foreach INTEGRAL,$(INTEGRALS),lib/lib$(INTEGRAL).so) : lib/lib%%.so : %%_weighted_integral.o
	$(MAKE) -C ../$* static
	$(CXX) -o $@ -shared $*_weighted_integral.o ../$*/lib$*.a $(LDFLAGS)

ifdef SECDEC_WITH_CUDA_FLAGS
lib/lib$(NAME).so : amplitude.o
	$(CXX) -o $@ -shared $+ $(LDFLAGS)
else
lib/lib$(NAME).so : amplitude.o
	$(CXX) -o $@ -shared $+ $(LDFLAGS) -Wl,-undefined,dynamic_lookup
endif

$(foreach INTEGRAL,$(INTEGRALS),../lib/lib$(INTEGRAL).so) ../lib/lib$(NAME).so : ../lib/%% : lib/%%
	mv $< $@

clean :
	rm -f *.o *.so *.a

very-clean : clean
