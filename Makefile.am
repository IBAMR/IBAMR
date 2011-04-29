## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = doc ibtk src lib bin examples
ACLOCAL_AMFLAGS = -I m4

## Standard make targets.
lib: all
examples: lib
	@(cd examples && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
install-exec-local:
	@echo "===================================================================="
	@echo "===================================================================="
	@echo "WARNING: make install is not officially supported by IBAMR!"
	@echo "         use at your own risk!"
	@echo "===================================================================="
	@echo "===================================================================="

## Add various "extra" files to the IBAMR distribution.
EXTRA_DIST  = $(top_srcdir)/COPYRIGHT
EXTRA_DIST += $(top_srcdir)/autogen.sh
EXTRA_DIST += $(top_srcdir)/doc/Doxyfile.in
EXTRA_DIST += $(top_srcdir)/doc/footer.html
EXTRA_DIST += $(top_srcdir)/doc/header.html
EXTRA_DIST += $(top_srcdir)/doc/libstdc++.tag
EXTRA_DIST += $(top_srcdir)/doc/mainpage.dox
EXTRA_DIST += $(top_srcdir)/doc/online_Doxyfile.in
EXTRA_DIST += $(top_srcdir)/doc/samrai.tag
EXTRA_DIST += $(top_srcdir)/doc/stylesheet.css
EXTRA_DIST += $(top_srcdir)/include/ibamr/AdvDiffHierarchyIntegrator.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/AdvDiffHypPatchOps.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/AdvectHypPatchOps.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/CCHierarchyProjector.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/GodunovAdvector.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/HierarchyProjector.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBAnchorPointSpec.I
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBAnchorPointSpec.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBAnchorPointSpecFactory.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBBeamForceGen.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBBeamForceSpec.I
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBBeamForceSpec.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBBeamForceSpecFactory.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBDataPostProcessor.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBEulerianForceFunction.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBEulerianSourceFunction.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBFEHierarchyIntegrator.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBHDF5Initializer.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBHierarchyIntegrator.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBImplicitHierarchyIntegrator.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBImplicitJacobian.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBImplicitModHelmholtzOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBImplicitModHelmholtzPETScLevelSolver.I
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBImplicitModHelmholtzPETScLevelSolver.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBImplicitOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBImplicitSFROperator.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBImplicitSJROperator.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBInstrumentPanel.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBInstrumentationSpec.I
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBInstrumentationSpec.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBInstrumentationSpecFactory.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBKirchhoffRodForceGen.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBLagrangianForceStrategy.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBLagrangianForceStrategySet.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBLagrangianSourceStrategy.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBMovingTargetPointForceGen.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBMovingTargetPointForceSpec.I
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBMovingTargetPointForceSpec.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBMovingTargetPointForceSpecFactory.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBRodForceSpec.I
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBRodForceSpec.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBRodForceSpecFactory.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBSourceSpec.I
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBSourceSpec.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBSourceSpecFactory.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBSpringForceGen.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBSpringForceSpec.I
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBSpringForceSpec.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBSpringForceSpecFactory.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBStaggeredHierarchyIntegrator.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBStandardForceGen.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBStandardSourceGen.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBStandardInitializer.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBTargetPointForceGen.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBTargetPointForceSpec.I
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBTargetPointForceSpec.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/IBTargetPointForceSpecFactory.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSCoefs.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSHierarchyIntegrator.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSIntermediateVelocityBcCoef.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSProjectionBcCoef.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSStaggeredBlockFactorizationPreconditioner.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSStaggeredBoxRelaxationFACOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSStaggeredCenteredConvectiveOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSStaggeredHierarchyIntegrator.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSStaggeredIntermediateVelocityBcCoef.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSStaggeredPPMConvectiveOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSStaggeredPhysicalBoundaryHelper.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSStaggeredPressureBcCoef.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSStaggeredProjectionBcCoef.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSStaggeredProjectionPreconditioner.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSStaggeredStokesOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSStaggeredVCStokesOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/INSStaggeredVelocityBcCoef.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/ibamr_enums.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/ibamr_utilities.h
EXTRA_DIST += $(top_srcdir)/include/ibamr/namespaces.h

## Automatically link header files in $(top_srcdir)/include and
## refresh emacs TAGS file.
##
## Note that these are for developer use only.
copy_headers:
	cd ibtk && $(MAKE) relink_headers
	cd $(top_srcdir) && src/tools/includes.pl
relink_headers:
	cd ibtk && $(MAKE) relink_headers
	cd $(top_srcdir) && src/tools/includes.pl --link
refresh_tags:
	cd ibtk && $(MAKE) refresh_tags
	cd $(top_srcdir) && find . -name "*.[ChIf]" -print | etags -
