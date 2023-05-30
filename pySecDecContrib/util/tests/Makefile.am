check_PROGRAMS = test_integrator test_cuba_integrators test_cquad test_qmc test_series test_integrand_container test_deep_apply test_uncertainties test_amplitude test_coefficient_parser

AM_CPPFLAGS = -I$(top_srcdir)
if SECDEC_WITH_CUDA
AM_CPPFLAGS += -x cu -DSECDEC_WITH_CUDA=1
else
AM_LDFLAGS = -pthread
endif

SECDEC_CONTRIB = `python3 -m pySecDecContrib --dirname`

test_integrator_SOURCES = catch_amalgamated.cpp test_integrator.cpp catch_amalgamated.hpp
test_cuba_integrators_SOURCES = catch_amalgamated.cpp test_cuba_integrators.cpp catch_amalgamated.hpp
test_cquad_SOURCES = catch_amalgamated.cpp test_cquad.cpp catch_amalgamated.hpp
test_qmc_SOURCES = catch_amalgamated.cpp test_qmc.cpp catch_amalgamated.hpp
test_series_SOURCES = catch_amalgamated.cpp test_series.cpp catch_amalgamated.hpp
test_integrand_container_SOURCES = catch_amalgamated.cpp test_integrand_container.cpp catch_amalgamated.hpp
test_deep_apply_SOURCES = catch_amalgamated.cpp test_deep_apply.cpp catch_amalgamated.hpp
test_uncertainties_SOURCES = catch_amalgamated.cpp test_uncertainties.cpp catch_amalgamated.hpp
test_amplitude_SOURCES = catch_amalgamated.cpp test_amplitude.cpp catch_amalgamated.hpp
test_coefficient_parser_SOURCES = catch_amalgamated.cpp test_coefficient_parser.cpp catch_amalgamated.hpp

test_cuba_integrators_LDADD = -lcuba
test_cuba_integrators_CXXFLAGS = -I$(SECDEC_CONTRIB)/include
test_cuba_integrators_LDADD += -L$(SECDEC_CONTRIB)/lib

test_qmc_LDADD = -lgsl -lgslcblas
test_qmc_LDFLAGS = -pthread
test_qmc_CXXFLAGS = -I$(SECDEC_CONTRIB)/include
test_qmc_LDADD += -L$(SECDEC_CONTRIB)/lib

test_cquad_LDADD = -lgsl -lgslcblas
test_cquad_CXXFLAGS = -I$(SECDEC_CONTRIB)/include
test_cquad_LDADD += -L$(SECDEC_CONTRIB)/lib

test_amplitude_LDADD = -lgsl -lgslcblas -lcuba
test_amplitude_CXXFLAGS = -I$(SECDEC_CONTRIB)/include
test_amplitude_LDADD += -L$(SECDEC_CONTRIB)/lib

test_coefficient_parser_LDADD = -lgmp
test_coefficient_parser_CXXFLAGS = -I$(SECDEC_CONTRIB)/include -DSECDEC_CONTRIB=\"$(SECDEC_CONTRIB)\"
test_coefficient_parser_LDADD += -L$(SECDEC_CONTRIB)/lib

TESTS = $(check_PROGRAMS)
