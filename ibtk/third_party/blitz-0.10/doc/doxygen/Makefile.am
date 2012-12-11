

# Support doxygen documentation ?

docdir = $(datadir)/doc/$(PACKAGE)-$(VERSION)/doxygen

DISTCLEANFILES = Doxyfile doxygen-warning 

if DOXYGEN_DOC


html/index.html: Doxyfile
	@DOXYGEN@ 

# Document Targets
html: html/index.html
dvi: $(PACKAGE).dvi
pdf: $(PACKAGE).pdf
ps: $(PACKAGE).ps

$(PACKAGE).dvi: html/index.html
	@if test -d ./latex; then \
	echo "Making DVI manual"; \
	if test -f $@; then \
		$(RM) $@; \
	fi; \
		$(MAKE) -C ./latex refman.dvi; \
		$(LN_S) ./latex/refman.dvi $@; \
	fi

$(PACKAGE).pdf: html/index.html
	@if test -d ./latex; then \
	echo "Making PDF manual"; \
	if test -f $@; then \
		$(RM) $@; \
	fi; \
		$(MAKE) -C ./latex refman.pdf; \
		$(LN_S) ./latex/refman.pdf $@; \
	fi

$(PACKAGE).ps: html/index.html
	@if test -d ./latex; then \
	echo "Making PS manual"; \
	if test -f $@; then \
		$(RM) $@; \
	fi; \
		$(MAKE) -C ./latex refman.ps; \
		$(LN_S) ./latex/refman.pdf $@; \
	fi

doc: html/index.html
	@if test -d latex; then \
		$(MAKE) ps; \
		$(MAKE) pdf; \
	fi

else
doc:
endif


all-local: doc


# documentation distribution
doc-dist: doc
	@if test -d html; then \
		echo "please wait while archiving the docs."; \
		tar cf - -C html . | bzip2 --best -c > $(PACKAGE)-docs.tar.bz2; \
	fi

# Clean hooks
clean-local:
	@if test -f $(PACKAGE).pdf; then \
		$(RM) -f $(PACKAGE).pdf; \
	fi
	@if test -f $(PACKAGE).ps; then \
		$(RM) -f $(PACKAGE).ps; \
	fi
	@if test -d html; then \
		$(RM) -rf html; \
	fi
	@if test -d latex; then \
		$(RM) -rf latex; \
	fi

# Install hooks
make-install-dirs:
	@if test '!' -d $(DESTDIR)$(docdir) ; \
		then mkdir -p $(DESTDIR)$(docdir) ; \
	fi

install-data-hook: make-install-dirs
	@if test -d html; then \
		echo Installing HTML documentation into $(DESTDIR)$(docdir) ; \
		$(INSTALL_DATA) html/* $(DESTDIR)$(docdir) ; \
	fi
	@if test -f $(PACKAGE).pdf; then \
		echo Installing PDF manual into $(DESTDIR)$(docdir) ; \
		$(INSTALL_DATA) $(PACKAGE).pdf $(DESTDIR)$(docdir) ; \
	fi
	@if test -f $(PACKAGE).ps; then \
		echo Installing PS manual into $(DESTDIR)$(docdir) ; \
		$(INSTALL_DATA) $(PACKAGE).ps $(DESTDIR)$(docdir) ; \
	fi

uninstall-hook: 
	@echo Uninstalling doxygen documentation directory ; \
	rm -rf "$(DESTDIR)$(docdir)";



