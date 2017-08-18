# Setup
# =====
.DEFAULT_GOAL=runtests

# Running Tests
# =============

# Running all
# -----------
ALL_EXAMPLES=$(shell ls -d *[!test_suite]*/ | sed s\#/\#\#)
.PHONY : runtests
runtests : $(ALL_EXAMPLES)

# Runing selection
# ----------------
SELECTED_TESTS = difference bubble1L_dot_rank4 userdefined_cpp
.PHONY : runselected
runselected : $(SELECTED_TESTS)

# General rule
# ------------
.PHONY .SILENT : $(ALL_EXAMPLES)
$(ALL_EXAMPLES) : test_suite/catch.o
	echo ------------------------------------------------------------
	echo "START `date +'%s'`"
	echo ------------------------------------------------------------
	echo "TEST $@"
	echo ------------------------------------------------------------
	printf "Running 'make very-clean' in $@ ...\n"
	$(MAKE) -C "$@" very-clean \
	|| printf "make very-clean failed\n@@@ FAILURE @@@\n" \
		> "$@/test/test.log" \
	printf "Running 'make test' in $@...\n"
	$(MAKE) -C "$@" test \
	|| printf "make test failed\n@@@ FAILURE @@@\n" \
		> "$@/test/test.log"
	echo ------------------------------------------------------------
	echo "DONE `date +'%s'`"
	echo ------------------------------------------------------------

test_suite/catch.o :
	$(MAKE) -C test_suite catch.o

# Printing Results
# ================
.SILENT .PHONY : summarize
summarize :
	printf "\n"
	printf "\033[0;32m"
	find . -type d -name 'test' -exec grep -Hl "@@@ SUCCESS @@@" \{}/test.log \
		   2>/dev/null \; | \
		sed 's/\.\/\([A-Za-z0-9_-]*\)\/test\/test.log/+ \1 (succeeded)/' | \
		sort
	printf "\033[0m"
	printf "\033[1;31m"
	find . -type d -name 'test' -exec grep -Hl "@@@ FAILURE @@@" \{}/test.log \
		   2>/dev/null \; | \
		sed 's/\.\/\([A-Za-z0-9_-]*\)\/test\/test.log/+ \1 (FAILED)/' | \
		sort
	find . -type d -name 'test' -exec grep -HL "@@@ FAILURE\|SUCCESS @@@" \{}/test.log \; | \
		sed 's/\.\/\([A-Za-z0-9_-]*\)\/test\/test.log/+ \1 (UNKNOWN)/' | \
		sort
	printf "\033[0m"
	printf "\n"

	# return nonzero exit status if any example does not report success
	for dir in $(ALL_EXAMPLES); do test -e $$dir/test/test.log || exit 1 && grep -HL "@@@ SUCCESS @@@" $$dir/test/test.log 2>&1 >/dev/null || exit 1 ; done

.SILENT .PHONY : summarize-selected
summarize-selected :
	printf "\n"
	printf "\033[0;32m"
	find $(SELECTED_TESTS) -type d -name 'test' -exec grep -Hl "@@@ SUCCESS @@@" \{}/test.log \
		   2>/dev/null \; | \
		sed 's/\([A-Za-z0-9_-]*\)\/test\/test.log/+ \1 (succeeded)/' | \
		sort
	printf "\033[0m"
	printf "\033[1;31m"
	find $(SELECTED_TESTS) -type d -name 'test' -exec grep -Hl "@@@ FAILURE @@@" \{}/test.log \
		   2>/dev/null \; | \
		sed 's/\([A-Za-z0-9_-]*\)\/test\/test.log/+ \1 (FAILED)/' | \
		sort
	find $(SELECTED_TESTS) -type d -name 'test' -exec grep -HL "@@@ FAILURE\|SUCCESS @@@" \{}/test.log \; | \
		sed 's/\([A-Za-z0-9_-]*\)\/test\/test.log/+ \1 (UNKNOWN)/' | \
		sort
	printf "\033[0m"
	printf "\n"

	# return nonzero exit status if any example does not report success
	for dir in $(SELECTED_TESTS); do test -e $$dir/test/test.log || exit 1 && grep -HL "@@@ SUCCESS @@@" $$dir/test/test.log 2>&1 >/dev/null || exit 1 ; done

# Cleaning up
# ===========
.PHONY : clean
clean :
	for f in * ; \
	do \
	if test -d "$$f" -a -f "$$f/Makefile" ; \
	then \
		$(MAKE) -C $$f very-clean ; \
	fi \
	done

very-clean : clean
