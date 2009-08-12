## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = doc ibtk src
ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST  = autogen.sh
EXTRA_DIST += doc/Doxyfile
EXTRA_DIST += doc/libstdc++.tag
EXTRA_DIST += include/ibamr/AdvDiffHierarchyIntegrator.h
EXTRA_DIST += include/ibamr/AdvDiffHypPatchOps.h
EXTRA_DIST += include/ibamr/AdvectHypPatchOps.h
EXTRA_DIST += include/ibamr/GodunovAdvector.h
EXTRA_DIST += include/ibamr/HierarchyProjector.h
EXTRA_DIST += include/ibamr/IBAnchorPointSpec.I
EXTRA_DIST += include/ibamr/IBAnchorPointSpec.h
EXTRA_DIST += include/ibamr/IBAnchorPointSpecFactory.h
EXTRA_DIST += include/ibamr/IBBeamForceGen.h
EXTRA_DIST += include/ibamr/IBBeamForceSpec.I
EXTRA_DIST += include/ibamr/IBBeamForceSpec.h
EXTRA_DIST += include/ibamr/IBBeamForceSpecFactory.h
EXTRA_DIST += include/ibamr/IBDataPostProcessor.h
EXTRA_DIST += include/ibamr/IBEulerianForceSetter.h
EXTRA_DIST += include/ibamr/IBEulerianSourceSetter.h
EXTRA_DIST += include/ibamr/IBHDF5Initializer.h
EXTRA_DIST += include/ibamr/IBHierarchyIntegrator.h
EXTRA_DIST += include/ibamr/IBImplicitHierarchyIntegrator.h
EXTRA_DIST += include/ibamr/IBInstrumentPanel.h
EXTRA_DIST += include/ibamr/IBInstrumentationSpec.I
EXTRA_DIST += include/ibamr/IBInstrumentationSpec.h
EXTRA_DIST += include/ibamr/IBInstrumentationSpecFactory.h
EXTRA_DIST += include/ibamr/IBLagrangianForceStrategy.h
EXTRA_DIST += include/ibamr/IBLagrangianForceStrategySet.h
EXTRA_DIST += include/ibamr/IBLagrangianSourceStrategy.h
EXTRA_DIST += include/ibamr/IBMovingTargetPointForceGen.h
EXTRA_DIST += include/ibamr/IBMovingTargetPointForceSpec.I
EXTRA_DIST += include/ibamr/IBMovingTargetPointForceSpec.h
EXTRA_DIST += include/ibamr/IBMovingTargetPointForceSpecFactory.h
EXTRA_DIST += include/ibamr/IBSpringForceGen.h
EXTRA_DIST += include/ibamr/IBSpringForceSpec.I
EXTRA_DIST += include/ibamr/IBSpringForceSpec.h
EXTRA_DIST += include/ibamr/IBSpringForceSpecFactory.h
EXTRA_DIST += include/ibamr/IBStaggeredHierarchyIntegrator.h
EXTRA_DIST += include/ibamr/IBStandardForceGen.h
EXTRA_DIST += include/ibamr/IBStandardInitializer.h
EXTRA_DIST += include/ibamr/IBTargetPointForceGen.h
EXTRA_DIST += include/ibamr/IBTargetPointForceSpec.I
EXTRA_DIST += include/ibamr/IBTargetPointForceSpec.h
EXTRA_DIST += include/ibamr/IBTargetPointForceSpecFactory.h
EXTRA_DIST += include/ibamr/INSCoefs.h
EXTRA_DIST += include/ibamr/INSHierarchyIntegrator.h
EXTRA_DIST += include/ibamr/INSIntermediateVelocityBcCoef.h
EXTRA_DIST += include/ibamr/INSProjectionBcCoef.h
EXTRA_DIST += include/ibamr/INSStaggeredBlockFactorizationPreconditioner.h
EXTRA_DIST += include/ibamr/INSStaggeredCenteredConvectiveOperator.h
EXTRA_DIST += include/ibamr/INSStaggeredHierarchyIntegrator.h
EXTRA_DIST += include/ibamr/INSStaggeredIntermediateVelocityBcCoef.h
EXTRA_DIST += include/ibamr/INSStaggeredPPMConvectiveOperator.h
EXTRA_DIST += include/ibamr/INSStaggeredPhysicalBoundaryHelper.h
EXTRA_DIST += include/ibamr/INSStaggeredPressureBcCoef.h
EXTRA_DIST += include/ibamr/INSStaggeredProjectionBcCoef.h
EXTRA_DIST += include/ibamr/INSStaggeredProjectionPreconditioner.h
EXTRA_DIST += include/ibamr/INSStaggeredStokesOperator.h
EXTRA_DIST += include/ibamr/INSStaggeredVelocityBcCoef.h
EXTRA_DIST += include/ibamr/TGACoefs.h
EXTRA_DIST += include/ibamr/get_version.h
EXTRA_DIST += sample/IB/input2d
EXTRA_DIST += sample/IB/input3d
EXTRA_DIST += sample/IB/main.C
EXTRA_DIST += sample/IB/Makefile
EXTRA_DIST += sample/IB/shell2d_64.spring
EXTRA_DIST += sample/IB/shell2d_64.vertex
EXTRA_DIST += sample/IB/sphere3d.spring
EXTRA_DIST += sample/IB/sphere3d.vertex

## Copy the sample sources to the build directory.
all-local:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(top_srcdir)/sample/IB/input[23]d $(top_builddir)/sample/IB ; \
	  cp -f $(top_srcdir)/sample/IB/main.C $(top_builddir)/sample/IB ; \
	  cp -f $(top_srcdir)/sample/IB/shell2d_64.{spring,vertex} $(top_builddir)/sample/IB ; \
	  cp -f $(top_srcdir)/sample/IB/sphere3d.{spring,vertex} $(top_builddir)/sample/IB ; \
	fi ;


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
