## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = stools doc src
ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST  = autogen.sh
EXTRA_DIST += bin/includes
EXTRA_DIST += doc/Doxyfile
EXTRA_DIST += doc/libstdc++.tag
EXTRA_DIST += include/ibamr/AdvDiffHierarchyIntegrator.h
EXTRA_DIST += include/ibamr/AdvDiffHypPatchOps.h
EXTRA_DIST += include/ibamr/ConvergenceMonitor.h
EXTRA_DIST += include/ibamr/get_version.h
EXTRA_DIST += include/ibamr/GodunovAdvector.h
EXTRA_DIST += include/ibamr/GodunovHypPatchOps.h
EXTRA_DIST += include/ibamr/HierarchyProjector.h
EXTRA_DIST += include/ibamr/IBEulerianForceSetter.h
EXTRA_DIST += include/ibamr/IBEulerianSourceSetter.h
EXTRA_DIST += include/ibamr/IBHierarchyIntegrator.h
EXTRA_DIST += include/ibamr/IBLagrangianForceStrategy.h
EXTRA_DIST += include/ibamr/IBLagrangianSourceStrategy.h
EXTRA_DIST += include/ibamr/IBLogicalCartMeshInitializer.h
EXTRA_DIST += include/ibamr/IBStandardInitializer.h
EXTRA_DIST += include/ibamr/INSHierarchyIntegrator.h
EXTRA_DIST += include/ibamr/LagSiloDataWriter.h
EXTRA_DIST += include/ibamr/LDataManager.h
EXTRA_DIST += include/ibamr/LEInteractor.h
EXTRA_DIST += include/ibamr/LNodeIndexDataFactory.h
EXTRA_DIST += include/ibamr/LNodeIndexData.h
EXTRA_DIST += include/ibamr/LNodeIndexData.I
EXTRA_DIST += include/ibamr/LNodeIndex.h
EXTRA_DIST += include/ibamr/LNodeIndex.I
EXTRA_DIST += include/ibamr/LNodeIndexSet.h
EXTRA_DIST += include/ibamr/LNodeIndexSet.I
EXTRA_DIST += include/ibamr/LNodeIndexVariable.h
EXTRA_DIST += include/ibamr/LNodeJacobianInitStrategy.h
EXTRA_DIST += include/ibamr/LNodeLevelData.h
EXTRA_DIST += include/ibamr/LNodeLevelData.I
EXTRA_DIST += include/ibamr/LNodePosnInitStrategy.h
EXTRA_DIST += include/ibamr/PhysicalBCDataStrategy.h
EXTRA_DIST += include/ibamr/SetDataStrategy.h
EXTRA_DIST += include/ibamr/SpringForceGen.h
EXTRA_DIST += include/ibamr/SpringForceSpecFactory.h
EXTRA_DIST += include/ibamr/SpringForceSpec.h
EXTRA_DIST += include/ibamr/StashableFactory.h
EXTRA_DIST += include/ibamr/Stashable.h
EXTRA_DIST += include/ibamr/StashableManager.h
EXTRA_DIST += include/ibamr/StashableManager.I
EXTRA_DIST += include/ibamr/StashableStream.h
EXTRA_DIST += include/ibamr/StashableStream.I
EXTRA_DIST += include/ibamr/TargetPointForceGen.h
EXTRA_DIST += include/ibamr/TargetPointForceSpecFactory.h
EXTRA_DIST += include/ibamr/TargetPointForceSpec.h
EXTRA_DIST += include/ibamr/TargetPointInitializer.h
EXTRA_DIST += sample/IB/input2d
EXTRA_DIST += sample/IB/input3d
EXTRA_DIST += sample/IB/main.C
EXTRA_DIST += sample/IB/shell2d_64.edges
EXTRA_DIST += sample/IB/shell2d_64.vertices
EXTRA_DIST += sample/IB/sphere3d.edges
EXTRA_DIST += sample/IB/sphere3d.vertices

## Copy the sample sources to the build directory.
all-local:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(top_srcdir)/sample/IB/input[23]d $(top_builddir)/sample/IB ; \
	  cp -f $(top_srcdir)/sample/IB/main.C $(top_builddir)/sample/IB ; \
	  cp -f $(top_srcdir)/sample/IB/shell2d_64.{edges,vertices} $(top_builddir)/sample/IB ; \
	  cp -f $(top_srcdir)/sample/IB/sphere3d.{edges,vertices} $(top_builddir)/sample/IB ; \
	fi ;


## Automatically link header files in $(top_srcdir)/include and
## refresh emacs TAGS file.
##
## Note that these are for developer use only.
copy_headers:
	cd stools && $(MAKE) relink_headers
	cd $(top_srcdir) && bin/includes
relink_headers:
	cd stools && $(MAKE) relink_headers
	cd $(top_srcdir) && bin/includes --link
refresh_tags:
	cd stools && $(MAKE) refresh_tags
	cd $(top_srcdir) && find . -name "*.[ChIf]" -print | etags -
