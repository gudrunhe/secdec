check_PROGRAMS = test_integrator test_cuba_integrators test_cquad test_series test_integrand_container test_deep_apply test_uncertainties

AM_CPPFLAGS = -I$(top_srcdir)
if SECDEC_WITH_CUDA
AM_CPPFLAGS += -DSECDEC_WITH_CUDA=1
endif

test_integrator_SOURCES = main.cpp test_integrator.cpp catch.hpp
test_cuba_integrators_SOURCES = main.cpp test_cuba_integrators.cpp catch.hpp
test_cquad_SOURCES = main.cpp test_cquad.cpp catch.hpp
test_series_SOURCES = main.cpp test_series.cpp catch.hpp
test_integrand_container_SOURCES = main.cpp test_integrand_container.cpp catch.hpp
test_deep_apply_SOURCES = main.cpp test_deep_apply.cpp catch.hpp
test_uncertainties_SOURCES = main.cpp test_uncertainties.cpp catch.hpp

test_cuba_integrators_LDADD = -lcuba
if HAVE_SECDEC_CONTRIB
test_cuba_integrators_CXXFLAGS = -I$(SECDEC_CONTRIB)/include
test_cuba_integrators_LDADD += -L$(SECDEC_CONTRIB)/lib
endif

test_cquad_LDADD = -lgsl -lgslcblas
if HAVE_SECDEC_CONTRIB
test_cquad_CXXFLAGS = -I$(SECDEC_CONTRIB)/include
test_cquad_LDADD += -L$(SECDEC_CONTRIB)/lib
endif

TESTS = $(check_PROGRAMS)
