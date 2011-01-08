LATEX = latex
DVIPS = dvips
DVIPDF = dvipdf
BIBTEX = bibtex

EPSFIGS = nine-two-sum.eps normal_sum_prod.eps qd_add.eps qd_add_proof.eps \
          qd_add_qd_d.eps qd_mul_accum.eps qd_mul_qd_d.eps quick-two-sum.eps \
          six-three-sum.eps three-sum-2.eps three-sum-3.eps three-sum.eps \
          two-prod.eps two-sum.eps 

all: qd.pdf

dvi: qd.dvi

pdf: qd.pdf

ps: qd.ps

qd.dvi: qd.tex qd.bib $(EPSFIGS)
	$(LATEX) $< && $(BIBTEX) qd && \
	$(LATEX) $< && $(LATEX) $<

qd.ps: qd.dvi
	$(DVIPS) -o $@ $<

qd.pdf: qd.dvi
	$(DVIPDF) $<

clean:
	rm -f qd.dvi qd.blg qd.bbl qd.toc qd.log qd.aux 

maintainer-clean: clean
	rm -f qd.ps qd.pdf

.PHONY: clean maintainer-clean all

