ACLOCAL_AMFLAGS = -I m4
SUBDIRS = config src include tests fortran

BUILT_SOURCES = docs/qd.pdf
dist_doc_DATA = README docs/qd.pdf
dist_noinst_DATA = docs/Makefile \
  docs/qd.tex docs/qd.bib \
  docs/nine-two-sum.eps docs/nine-two-sum.fig \
  docs/normal_sum_prod.eps docs/normal_sum_prod.fig \
  docs/qd_add.eps docs/qd_add.fig \
  docs/qd_add_proof.eps docs/qd_add_proof.fig \
  docs/qd_add_qd_d.eps docs/qd_add_qd_d.fig \
  docs/qd_mul_accum.eps docs/qd_mul_accum.fig \
  docs/qd_mul_qd_d.eps docs/qd_mul_qd_d.fig \
  docs/quick-two-sum.eps docs/quick-two-sum.fig \
  docs/six-three-sum.eps docs/six-three-sum.fig \
  docs/three-sum-2.eps docs/three-sum-2.fig \
  docs/three-sum-3.eps docs/three-sum-3.fig \
  docs/three-sum.eps docs/three-sum.fig \
  docs/two-prod.eps docs/two-prod.fig \
  docs/two-sum.eps docs/two-sum.fig

doc: docs/qd.pdf

docs/qd.pdf:
	$(MAKE) -C docsqd.pdf

changelog:
	git log >ChangeLog

dist-hook:
	rm -f ${distdir}/ChangeLog && git log >${distdir}/ChangeLog

cpp-demo:
	$(MAKE) -C tests demo

if HAVE_FORTRAN

fortran-demo:
	$(MAKE) -C fortran demo

demo: cpp-demo fortran-demo

else

fortran-demo:
	@echo "You need a Fortran 95 compiler to build fortran demo programs." && false

demo: cpp-demo

endif

time:
	$(MAKE) -C tests time

bin_SCRIPTS=qd-config

.PHONY: changelog cpp-demo fortran-demo demo time doc
