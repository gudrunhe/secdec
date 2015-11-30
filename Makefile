# This Makefile implements common tasks needed by developers
# A list of implemented rules can be obtained by the command "make help"

# Set the timeout for a test case to finish here
TIMEOUT=120


.DEFAULT_GOAL=check
.PHONY .SILENT : help
help :
	echo
	echo "    Implemented targets:"
	echo
	echo "    check        use nosetests to test with python 2.7 and 3"
	echo "    checkX       use nosetests to test with python 2.7 or 3,"
	echo "                 where X is one of {2,3}"
	echo "    active_check use nosetests to run only tests marked as"
	echo "                 \"active\" using nosetests-2.7 and nosetests3"
	echo "    fast_check   use nosetests to run only quick tests"
	echo "                 using nosetests-2.7 and nosetests3"
	echo "    clean        delete compiled and temporary files"
	echo "    coverage     produce and show a code coverage report"
	echo "    doc          build the html documentation using sphinx"
	echo "    doc-pdf      build the pdf documentation using sphinx"
	echo "    help         show this message"
	echo "    run-examples run all examples using python 2 and 3"
	echo "    show-todos   show todo marks in the source code"
	echo

.PHONY : clean
clean:
	#remove build doc
	make -C ./doc clean

	#remove .pyc files created by python 2.7
	rm -f ./*.pyc
	find -P . -name '*.pyc' -delete

	#remove .pyc files crated by python 3
	rm -rf ./__pycache__
	find -P . -name __pycache__ -delete

	#remove backup files
	find -P . -name '*~' -delete

	#remove files created by coverage
	rm -f .coverage
	rm -rf coverage

.PHONY : check
check : check2 check3

.PHONY : check2
check2 :
	@ # run tests
	nosetests-2.7 --processes=-1 --process-timeout=$(TIMEOUT)

.PHONY : check3
check3 :
	@ # run tests
	nosetests3 --processes=-1 --process-timeout=$(TIMEOUT)

.PHONY : active_check
active_check :
	nosetests-2.7 -a 'active' --processes=-1 --process-timeout=$(TIMEOUT)
	nosetests3    -a 'active' --processes=-1 --process-timeout=$(TIMEOUT)

.PHONY : fast_check
fast_check :
	nosetests-2.7 -a '!slow' --processes=-1 --process-timeout=$(TIMEOUT)
	nosetests3    -a '!slow' --processes=-1 --process-timeout=$(TIMEOUT)

.PHONY : doc
doc :
	cd doc && make html

.PHONY : doc-pdf
doc-pdf :
	cd doc; make latexpdf

.PHONY : run-examples
run-examples :
	cd examples ; \
	for file in $$(ls) ; do \
	    echo running $${file} with python2 && \
	    python2 $${file} || exit 1 && \
	    echo running $${file} with python3 && \
	    python3 $${file} || exit 1 && \
	; \
	done

.SILENT .PHONY : show-todos
grep_cmd  = grep -riG [^"au""sphinx.ext."]todo --color=auto --exclude=Makefile --exclude-dir=.git
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
	echo -e $(begin_red)"* o .git directories                                              *"$(end_red) ; \
	echo -e $(begin_red)"*******************************************************************"$(end_red) ; \
	echo

.PHONY : coverage
coverage :
	rm -rf coverage
	nosetests --with-coverage --cover-package=pySecDec --cover-html --cover-html-dir=coverage
	xdg-open coverage/index.html
