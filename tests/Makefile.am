include $(top_srcdir)/config/Make-rules
TEST_DIRECTORIES = adv_diff advect IB IBTK IBFE multiphase_flow physical_boundary refine spread wave_tank

tests: $(TEST_DIRECTORIES)

$(TEST_DIRECTORIES):
	cd $@ && $(MAKE) $(AM_MAKEFLAGS) tests
.PHONY: tests $(TEST_DIRECTORIES)
