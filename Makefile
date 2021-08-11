# This Makefile implements common tasks needed by developers
# A list of implemented rules can be obtained by the command "make help"

# timeout (in seconds) for all test cases
TIMEOUT=600

# Auto detect nvcc
ifeq ($(CXX), nvcc)
SECDEC_WITH_CUDA = 1
endif


.DEFAULT_GOAL=build
.PHONY .SILENT : help
help :
	echo
	echo "    Implemented targets:"
	echo
	echo "    dependencies install the Python dependencies (via pip)"
	echo "    build        build the contributed software and any needed"
	echo "                 autogenerated files (except docs)"
	echo "    check        run tests with Nose, run doctest using Sphinx,"
	echo "                 run the tests of the  SecDecUtil package,"
	echo "                 and run all examples"
	echo "    active-check use Nose to run only tests marked as \"active\""
	echo "    fast-check   use Nose to run only quick tests"
	echo "    util-check   run the tests of the SecDecUtil package"
	echo "    doctest      run doctest using Sphinx"
	echo "    dist         create distribution files for PyPI"
	echo "    clean        delete compiled and temporary files"
	echo "    coverage     produce and show a code coverage report"
	echo "    doc          run \"doc-html\" and \"doc-pdf\""
	echo "    doc-html     build the html documentation using sphinx"
	echo "    doc-pdf      build the pdf documentation using sphinx"
	echo "    help         show this message"
	echo "    high-level   run the high level tests"
	echo "    show-todos   show todo marks in the source code"
	echo

.PHONY : clean
clean:
	# remove build doc
	$(MAKE) -C ./doc clean

	# clean high level tests
	$(MAKE) -C ./high_level_tests clean

	# remove cpp doctests
	$(MAKE) -C ./doc/source/cpp_doctest clean

	# clean util
	if [ -f util/Makefile ] ; then $(MAKE) -C ./util clean ; fi

	# remove .pyc files created by python 3
	find -P . -name '*.pyc' -delete
	find -P . -name __pycache__ -delete

	# remove backup files
	find -P . -name '*~' -delete

	# remove Mac .DS_store files
	find -P . -name '.DS_Store' -delete

	# remove files created by coverage
	rm -f .coverage
	rm -rf coverage

	# remove egg-info
	rm -rf pySecDec.egg-info

	# remove build/ und dist/
	rm -rf build/ dist/

	# remove the SecDecUtil tarball
	rm -f util/secdecutil-*.tar.gz

	python3 -m SCons -cc

.PHONY: dependencies
dependencies:
	python3 -m pip install --user toml
	python3 -m pip install --user $$(python3 -c 'import toml; m = toml.load("pyproject.toml"); print(" ".join(m["build-system"]["requires"] + m["project"]["dependencies"] + sum(m["project"]["optional-dependencies"].values(), [])));')

.PHONY: build
build:
	python3 -m SCons build -j2

.PHONY : check
check : check3 doctest util-check

.PHONY : check3
check3 :
	@ # run tests
	python3 -m nose --processes=-1 --process-timeout=$(TIMEOUT)

.PHONY : active-check
active-check :
	python3 -m nose -a 'active' --processes=-1 --process-timeout=$(TIMEOUT)

.PHONY : fast-check
fast-check :
	python3 -m nose -a '!slow' --processes=-1 --process-timeout=$(TIMEOUT)

.PHONY : util-check
util-check :
	export SECDEC_WITH_CUDA=$(SECDEC_WITH_CUDA) && \
	export SECDEC_CONTRIB=$$(python3 -m pySecDecContrib --dirname) && \
	cd pySecDecContrib/util && \
	if [ -f Makefile ] ; \
	then \
		$(MAKE) check ; \
	else \
		autoreconf -i && \
		./configure --prefix=`pwd` && \
		$(MAKE) check ; \
	fi

.PHONY : doc
doc : doc-html doc-pdf

.PHONY : doc-html
doc-html :
	$(MAKE) -C doc html

.PHONY : doc-pdf
doc-pdf :
	$(MAKE) -C doc latexpdf

.PHONY : doctest
doctest :
	$(MAKE) -C doc doctest

.PHONY : high-level
high-level :
	# '$(MAKE) -C high_level_tests summarize' forwards the error if an example fails
	# 'exit 1' forwards the error out of the shell's for loop
	export DIRNAME=high_level_tests ; \
	for PYTHON in python3 ; do \
		PYTHON=$$PYTHON $(MAKE) -C $$DIRNAME && $(MAKE) -C $$DIRNAME summarize || exit 1 \
	; \
	done

.PHONY : dist
dist : build
	python3 -m SCons dist

.SILENT .PHONY : show-todos
grep_cmd  = grep -riG [^"au""sphinx.ext."]todo --color=auto --exclude=Makefile --exclude-dir=.git --exclude=catch.hpp
begin_red = "\033[0;31m"
end_red   = "\033[0m"
show-todos :
	# suppress errors here
	# note that no todo found is considered as error
	$(grep_cmd) . ; \
	echo ; 	echo ; \
	echo -e $(begin_red)"*******************************************************************"$(end_red) ; \
	echo -e $(begin_red)"* The following files and directories are NOT searched for TODOs: *"$(end_red) ; \
	echo -e $(begin_red)"* o makefiles                                                     *"$(end_red) ; \
	echo -e $(begin_red)"* o files named 'catch.hpp'                                       *"$(end_red) ; \
	echo -e $(begin_red)"* o .git directories                                              *"$(end_red) ; \
	echo -e $(begin_red)"*******************************************************************"$(end_red) ; \
	echo

.PHONY : coverage
coverage :
	rm -rf coverage
	python3 -m nose --with-coverage --cover-package=pySecDec --cover-html --cover-html-dir=coverage
	xdg-open coverage/index.html
