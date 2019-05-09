include $(top_srcdir)/config/Make-rules
TEST_DIRECTORIES = advect IBTK IBFE multiphase_flow

tests: $(TEST_DIRECTORIES)

$(TEST_DIRECTORIES):
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) tests
.PHONY: tests $(TEST_DIRECTORIES)
