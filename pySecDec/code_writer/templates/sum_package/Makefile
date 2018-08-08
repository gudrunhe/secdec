include Makefile.conf

# TODO : move below source after .DEFAULT_GOAL in Makefile.conf is fixed
integrate_$(NAME) : integrate_$(NAME).o $(foreach INTEGRAL,$(INTEGRALS),$(INTEGRAL)/lib$(INTEGRAL).a)
	$(CXX) $(CXXFLAGS) -o $@ $< $(patsubst %%,-L%%,$(INTEGRALS)) $(patsubst %%,-l%%,$(INTEGRALS)) $(LDFLAGS)

source :
	for INTEGRAL in $(INTEGRALS); do $(MAKE) -C $${INTEGRAL} $@; done

# rule to build the integral libraries
$(foreach INTEGRAL,$(INTEGRALS),$(INTEGRAL)/lib$(INTEGRAL).a) :
	$(MAKE) -C $(dir $@) $(notdir $@)

very-clean : clean
	for dir in */; do $(MAKE) -C $$dir $@; done

clean:
	for dir in */; do $(MAKE) -C $$dir $@; done
	rm -f *.o *.so *.a integrate_$(NAME)
