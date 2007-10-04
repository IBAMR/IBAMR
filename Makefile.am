## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = stools doc src
ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST  = autogen.sh
EXTRA_DIST += doc/Doxyfile
EXTRA_DIST += doc/libstdc++.tag
EXTRA_DIST += include/ibamr/AdvDiffHierarchyIntegrator.h
EXTRA_DIST += include/ibamr/AdvDiffHypPatchOps.h
EXTRA_DIST += include/ibamr/AdvectHypPatchOps.h
EXTRA_DIST += include/ibamr/get_version.h
EXTRA_DIST += include/ibamr/GodunovAdvector.h
EXTRA_DIST += include/ibamr/HierarchyProjector.h
EXTRA_DIST += include/ibamr/IBBeamForceGen.h
EXTRA_DIST += include/ibamr/IBBeamForceSpecFactory.h
EXTRA_DIST += include/ibamr/IBBeamForceSpec.h
EXTRA_DIST += include/ibamr/IBBeamForceSpec.I
EXTRA_DIST += include/ibamr/IBEulerianForceSetter.h
EXTRA_DIST += include/ibamr/IBEulerianSourceSetter.h
EXTRA_DIST += include/ibamr/IBHierarchyIntegrator.h
EXTRA_DIST += include/ibamr/IBImplicitHierarchyIntegrator.h
EXTRA_DIST += include/ibamr/IBInstrumentPanel.h
EXTRA_DIST += include/ibamr/IBInstrumentationSpecFactory.h
EXTRA_DIST += include/ibamr/IBInstrumentationSpec.h
EXTRA_DIST += include/ibamr/IBInstrumentationSpec.I
EXTRA_DIST += include/ibamr/IBLagrangianForceStrategy.h
EXTRA_DIST += include/ibamr/IBLagrangianForceStrategySet.h
EXTRA_DIST += include/ibamr/IBLagrangianSourceStrategy.h
EXTRA_DIST += include/ibamr/IBSpringForceGen.h
EXTRA_DIST += include/ibamr/IBSpringForceSpecFactory.h
EXTRA_DIST += include/ibamr/IBSpringForceSpec.h
EXTRA_DIST += include/ibamr/IBSpringForceSpec.I
EXTRA_DIST += include/ibamr/IBStandardForceGen.h
EXTRA_DIST += include/ibamr/IBStandardInitializer.h
EXTRA_DIST += include/ibamr/IBTargetPointForceGen.h
EXTRA_DIST += include/ibamr/IBTargetPointForceSpecFactory.h
EXTRA_DIST += include/ibamr/IBTargetPointForceSpec.h
EXTRA_DIST += include/ibamr/IBTargetPointForceSpec.I
EXTRA_DIST += include/ibamr/INSHierarchyIntegrator.h
EXTRA_DIST += include/ibamr/LagSiloDataWriter.h
EXTRA_DIST += include/ibamr/LDataManager.h
EXTRA_DIST += include/ibamr/LEInteractor.h
EXTRA_DIST += include/ibamr/LNodeIndexDataFactory2.h
EXTRA_DIST += include/ibamr/LNodeIndexData2.h
EXTRA_DIST += include/ibamr/LNodeIndexData2.I
EXTRA_DIST += include/ibamr/LNodeIndex.h
EXTRA_DIST += include/ibamr/LNodeIndex.I
EXTRA_DIST += include/ibamr/LNodeIndexSet.h
EXTRA_DIST += include/ibamr/LNodeIndexSet.I
EXTRA_DIST += include/ibamr/LNodeIndexVariable2.h
EXTRA_DIST += include/ibamr/LNodeInitStrategy.h
EXTRA_DIST += include/ibamr/LNodeLevelData.h
EXTRA_DIST += include/ibamr/LNodeLevelData.I
EXTRA_DIST += include/ibamr/NormalVelocityBcCoefAdapter.h
EXTRA_DIST += include/ibamr/ParabolicLevelIntegrator.h
EXTRA_DIST += include/ibamr/StashableFactory.h
EXTRA_DIST += include/ibamr/Stashable.h
EXTRA_DIST += include/ibamr/StashableManager.h
EXTRA_DIST += include/ibamr/StashableManager.I
EXTRA_DIST += include/ibamr/StashableStream.h
EXTRA_DIST += include/ibamr/StashableStream.I
EXTRA_DIST += sample/IB/input2d
EXTRA_DIST += sample/IB/input3d
EXTRA_DIST += sample/IB/main.C
EXTRA_DIST += sample/IB/shell2d_64.spring
EXTRA_DIST += sample/IB/shell2d_64.vertex
EXTRA_DIST += sample/IB/sphere3d.spring
EXTRA_DIST += sample/IB/sphere3d.vertex
EXTRA_DIST += src/tools/autostarter.pl
EXTRA_DIST += src/tools/farchiver.pl
EXTRA_DIST += src/tools/includes.pl
EXTRA_DIST += src/tools/scale_spring_rest_length.pl
EXTRA_DIST += src/tools/scale_spring_stiffness.pl
EXTRA_DIST += src/tools/set_spring_rest_length.pl
EXTRA_DIST += src/tools/set_spring_stiffness.pl

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
	cd stools && $(MAKE) relink_headers
	cd $(top_srcdir) && src/tools/includes.pl
relink_headers:
	cd stools && $(MAKE) relink_headers
	cd $(top_srcdir) && src/tools/includes.pl --link
refresh_tags:
	cd stools && $(MAKE) refresh_tags
	cd $(top_srcdir) && find . -name "*.[ChIf]" -print | etags -
