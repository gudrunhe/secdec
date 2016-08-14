check_PROGRAMS = test_series test_integrand_container test_algorithm test_uncertainties

AM_CPPFLAGS = -I$(top_srcdir)

test_series_SOURCES = main.cpp test_series.cpp catch.hpp
test_integrand_container_SOURCES = main.cpp test_integrand_container.cpp catch.hpp
test_algorithm_SOURCES = main.cpp test_algorithm.cpp catch.hpp
test_uncertainties_SOURCES = main.cpp test_uncertainties.cpp catch.hpp

TESTS = $(check_PROGRAMS)
