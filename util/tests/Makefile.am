check_PROGRAMS = test_integrators test_series test_integrand_container test_deep_apply test_uncertainties

AM_CPPFLAGS = -I$(top_srcdir)

test_integrators_SOURCES = main.cpp test_integrators.cpp catch.hpp
test_series_SOURCES = main.cpp test_series.cpp catch.hpp
test_integrand_container_SOURCES = main.cpp test_integrand_container.cpp catch.hpp
test_deep_apply_SOURCES = main.cpp test_deep_apply.cpp catch.hpp
test_uncertainties_SOURCES = main.cpp test_uncertainties.cpp catch.hpp

test_integrators_LDADD = -lcuba
if HAVE_SECDEC_CONTRIB
test_integrators_CXXFLAGS = -I$(SECDEC_CONTRIB)/include
test_integrators_LDADD += -L$(SECDEC_CONTRIB)/lib
endif

TESTS = $(check_PROGRAMS)
