.SECONDEXPANSION:  # to expand $$(@D)/.DIR
.SUFFIXES:         # Clear .SUFFIXES because we don't use implicit rules
.DELETE_ON_ERROR:  # Delete likely-corrupt target file if rule fails

OBJDIR := $(abspath obj)
LIBDIR := $(abspath lib)
BINDIR ?= bin
INCDIR ?= include

##### Workarounds for Broken Architectures #####
# old cygwin versions
ifeq ($(PETSC_CYGWIN_BROKEN_PIPE),1)
ifeq ($(shell basename $(AR)),ar)
  V ?=1
endif
endif
ifeq ($(V),)
  quiet_HELP := "Use \"$(MAKE) V=1\" to see the verbose compile lines.\n"
  quiet = @printf $(quiet_HELP)$(eval quiet_HELP:=)"  %10s %s\n" "$1$2" "$@"; $($1)
else ifeq ($(V),0)		# Same, but do not print any help
  quiet = @printf "  %10s %s\n" "$1$2" "$@"; $($1)
else				# Show the full command line
  quiet = $($1)
endif

##### Versioning #####
IBAMR_VERSION_MAJOR := $(shell awk '/\#define IBAMR_VERSION_MAJOR/{print $$3;}' $(IBAMR_DIR)/include/ibamrversion.h)
IBAMR_VERSION_MINOR := $(shell awk '/\#define IBAMR_VERSION_MINOR/{print $$3;}' $(IBAMR_DIR)/include/ibamrversion.h)
IBAMR_VERSION_SUBMINOR := $(shell awk '/\#define IBAMR_VERSION_SUBMINOR/{print $$3;}' $(IBAMR_DIR)/include/ibamrversion.h)
IBAMR_VERSION_RELEASE := $(shell awk '/\#define IBAMR_VERSION_RELEASE/{print $$3;}' $(IBAMR_DIR)/include/ibamrversion.h)
ifeq ($(IBAMR_VERSION_RELEASE),0)
  IBAMR_VERSION_MINOR := 0$(IBAMR_VERSION_MINOR)
endif
libibamr_abi_version := $(IBAMR_VERSION_MAJOR).$(IBAMR_VERSION_MINOR)
libibamr_lib_version := $(libibamr_abi_version).$(IBAMR_VERSION_SUBMINOR)

##### Functions #####
# Function to name shared library $(call SONAME_FUNCTION,libfoo,abiversion)
SONAME_FUNCTION ?= $(1).$(SL_LINKER_SUFFIX).$(2)
soname_function  = $(call SONAME_FUNCTION,$(1),$(libibamr_abi_version))
libname_function = $(call SONAME_FUNCTION,$(1),$(libibamr_lib_version))
# Function to link shared library $(call SL_LINKER_FUNCTION,libfoo,abiversion,libversion)
SL_LINKER_FUNCTION ?= -shared -Wl,-soname,$(call SONAME_FUNCTION,$(notdir $(1)),$(2))
basename_all   = $(basename $(basename $(basename $(basename $(1)))))
sl_linker_args = $(call SL_LINKER_FUNCTION,$(call basename_all,$@),$(libibamr_abi_version),$(libibamr_lib_version))
# Function to prefix directory that contains most recently-parsed makefile (current) if that directory is not ./
thisdir = $(addprefix $(dir $(lastword $(MAKEFILE_LIST))),$(1))
# Function to include makefile from subdirectories
incsubdirs = $(addsuffix /local.mk,$(call thisdir,$(1)))
# Function to change source filenames to object filenames
srctoobj = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(filter-out $(OBJDIR)/%,$(1)))

##### Libraries #####
libibamr_shared  := $(LIBDIR)/libibamr.$(SL_LINKER_SUFFIX)
libibamr_soname  := $(call soname_function,$(LIBDIR)/libibamr)
libibamr_libname := $(call libname_function,$(LIBDIR)/libibamr)
libibamr_static  := $(LIBDIR)/libibamr.$(AR_LIB_SUFFIX)
libibamr         := $(if $(filter-out no,$(BUILDSHAREDLIB)),$(libibamr_shared) $(libibamr_soname),$(libibamr_static))

##### Must define these, or thisdir does not work ######
srcs-core.cpp :=
srcs-ibtk.cpp :=
srcs-core.F   :=
srcs-ibtk.F   :=

srcs-IB_explicit_0.cpp :=
srcs-IB_explicit_1.cpp :=
srcs-IB_explicit_2.cpp :=
srcs-IB_explicit_3.cpp :=
srcs-IB_explicit_4.cpp :=
srcs-IB_explicit_5.cpp :=
srcs-IB_explicit_6.cpp :=
srcs-IB_implicit_0.cpp :=
srcs-IB_implicit_1.cpp :=
srcs-IB_implicit_2.cpp :=
srcs-IBFE_explicit_0.cpp :=
srcs-IBFE_explicit_1.cpp :=
srcs-IBFE_explicit_2.cpp :=
srcs-IBFE_explicit_3.cpp :=
srcs-IBFE_explicit_4.cpp :=
srcs-IBFE_explicit_5.cpp :=
srcs-IBFE_explicit_6.cpp :=
srcs-IBFE_explicit_7.cpp :=
srcs-IBFE_explicit_8.cpp :=
srcs-IMP_explicit_0.cpp :=
srcs-adv_diff_0.cpp :=
srcs-adv_diff_1.cpp :=
srcs-adv_diff_2.cpp :=
srcs-advect.cpp :=
srcs-navier_stokes_0.cpp :=
srcs-navier_stokes_1.cpp :=
srcs-navier_stokes_2.cpp :=
srcs-navier_stokes_3.cpp :=
srcs-navier_stokes_4.cpp :=
srcs-navier_stokes_5.cpp :=
srcs-navier_stokes_6.cpp :=
srcs-ConstraintIB_eel2d.cpp :=

#  IBFE_explicit_0 IBFE_explicit_1 IBFE_explicit_2 IBFE_explicit_3 IBFE_explicit_4 IBFE_explicit_5 IBFE_explicit_6 IBFE_explicit_7 IBFE_explicit_8
#  IMP_explicit_0
#  ConstraintIB_eel2d
examples := IB_explicit_0 IB_explicit_1 IB_explicit_2 IB_explicit_3 IB_explicit_4 IB_explicit_5 IB_explicit_6 IB_implicit_0 IB_implicit_1 IB_implicit_2 \
  adv_diff_0 adv_diff_1 adv_diff_2 advect \
  navier_stokes_0 advect navier_stokes_1 advect navier_stokes_2 advect navier_stokes_3 advect navier_stokes_4 advect navier_stokes_5 navier_stokes_6 \

##### Top level Rule #####
all : $(libibamr) $(examples)

##### Inclusions #####
# Recursively include files for all targets, needs to be defined before the source rules
include $(IBAMR_DIR)/local.mk

##### Rules #####
# make print VAR=the-variable
print:
	@echo $($(VAR))

IBAMR_COMPILE.c   = $(call quiet,$(cc_name)) -c $(PCC_FLAGS) $(CFLAGS) $(CCPPFLAGS) $(C_DEPFLAGS)
IBAMR_COMPILE.cxx = $(call quiet,CXX) -c $(CXX_FLAGS) $(CFLAGS) $(CCPPFLAGS) $(CXX_DEPFLAGS)
IBAMR_COMPILE.cu  = $(call quiet,CUDAC) -c $(CUDAC_FLAGS) --compiler-options="$(PCC_FLAGS) $(CXXFLAGS) $(CCPPFLAGS)"
IBAMR_GENDEPS.cu  = $(call quiet,CUDAC,.dep) --generate-dependencies $(CUDAC_FLAGS) --compiler-options="$(PCC_FLAGS) $(CXXFLAGS) $(CCPPFLAGS)"
IBAMR_COMPILE.F   = $(call quiet,FC) -c $(FC_FLAGS) $(FFLAGS) $(FCPPFLAGS) $(FC_DEPFLAGS)

pkgs := ibtk core
langs := c cu cpp F
concatlang = $(foreach lang, $(langs), $(srcs-$(1).$(lang):%.$(lang)=$(OBJDIR)/%.o))
srcs.o := $(foreach pkg, $(pkgs), $(call concatlang,$(pkg)))

srcs-IB_explicit_0.o := $(foreach ex, IB_explicit_0, $(call concatlang,$(ex)))
srcs-IB_explicit_1.o := $(foreach ex, IB_explicit_1, $(call concatlang,$(ex)))
srcs-IB_explicit_2.o := $(foreach ex, IB_explicit_2, $(call concatlang,$(ex)))
srcs-IB_explicit_3.o := $(foreach ex, IB_explicit_3, $(call concatlang,$(ex)))
srcs-IB_explicit_4.o := $(foreach ex, IB_explicit_4, $(call concatlang,$(ex)))
srcs-IB_explicit_5.o := $(foreach ex, IB_explicit_5, $(call concatlang,$(ex)))
srcs-IB_explicit_6.o := $(foreach ex, IB_explicit_6, $(call concatlang,$(ex)))
srcs-IB_implicit_0.o := $(foreach ex, IB_implicit_0, $(call concatlang,$(ex)))
srcs-IB_implicit_1.o := $(foreach ex, IB_implicit_1, $(call concatlang,$(ex)))
srcs-IB_implicit_2.o := $(foreach ex, IB_implicit_2, $(call concatlang,$(ex)))
srcs-IBFE_explicit_0.o := $(foreach ex, IBFE_explicit_0, $(call concatlang,$(ex)))
srcs-IBFE_explicit_1.o := $(foreach ex, IBFE_explicit_1, $(call concatlang,$(ex)))
srcs-IBFE_explicit_2.o := $(foreach ex, IBFE_explicit_2, $(call concatlang,$(ex)))
srcs-IBFE_explicit_3.o := $(foreach ex, IBFE_explicit_3, $(call concatlang,$(ex)))
srcs-IBFE_explicit_4.o := $(foreach ex, IBFE_explicit_4, $(call concatlang,$(ex)))
srcs-IBFE_explicit_5.o := $(foreach ex, IBFE_explicit_5, $(call concatlang,$(ex)))
srcs-IBFE_explicit_6.o := $(foreach ex, IBFE_explicit_6, $(call concatlang,$(ex)))
srcs-IBFE_explicit_7.o := $(foreach ex, IBFE_explicit_7, $(call concatlang,$(ex)))
srcs-IBFE_explicit_8.o := $(foreach ex, IBFE_explicit_8, $(call concatlang,$(ex)))
srcs-IMP_explicit_0.o := $(foreach ex, IMP_explicit_0, $(call concatlang,$(ex)))
srcs-adv_diff_0.o := $(foreach ex, adv_diff_0, $(call concatlang,$(ex)))
srcs-adv_diff_1.o := $(foreach ex, adv_diff_1, $(call concatlang,$(ex)))
srcs-adv_diff_2.o := $(foreach ex, adv_diff_2, $(call concatlang,$(ex)))
srcs-advect.o := $(foreach ex, advvect, $(call concatlang,$(ex)))
srcs-navier_stokes_0.o := $(foreach ex, navier_stokes_0, $(call concatlang,$(ex)))
srcs-navier_stokes_1.o := $(foreach ex, navier_stokes_1, $(call concatlang,$(ex)))
srcs-navier_stokes_2.o := $(foreach ex, navier_stokes_2, $(call concatlang,$(ex)))
srcs-navier_stokes_3.o := $(foreach ex, navier_stokes_3, $(call concatlang,$(ex)))
srcs-navier_stokes_4.o := $(foreach ex, navier_stokes_4, $(call concatlang,$(ex)))
srcs-navier_stokes_5.o := $(foreach ex, navier_stokes_5, $(call concatlang,$(ex)))
srcs-navier_stokes_6.o := $(foreach ex, navier_stokes_6, $(call concatlang,$(ex)))
srcs-ConstraintIB_eel2d.o := $(foreach ex, ConstraintIB_eel2d, $(call concatlang,$(ex)))

IB_explicit_0: $(srcs-IB_explicit_0.o)
	$(CLINKER) -o $@ $(srcs-IB_explicit_0.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IB_explicit_1: $(srcs-IB_explicit_1.o)
	$(CLINKER) -o $@ $(srcs-IB_explicit_1.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IB_explicit_2: $(srcs-IB_explicit_2.o)
	$(CLINKER) -o $@ $(srcs-IB_explicit_2.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IB_explicit_3: $(srcs-IB_explicit_3.o)
	$(CLINKER) -o $@ $(srcs-IB_explicit_3.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IB_explicit_4: $(srcs-IB_explicit_4.o)
	$(CLINKER) -o $@ $(srcs-IB_explicit_4.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IB_explicit_5: $(srcs-IB_explicit_5.o)
	$(CLINKER) -o $@ $(srcs-IB_explicit_5.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IB_explicit_6: $(srcs-IB_explicit_6.o)
	$(CLINKER) -o $@ $(srcs-IB_explicit_6.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IB_implicit_0: $(srcs-IB_implicit_0.o)
	$(CLINKER) -o $@ $(srcs-IB_implicit_0.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IB_implicit_1: $(srcs-IB_implicit_1.o)
	$(CLINKER) -o $@ $(srcs-IB_implicit_1.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IB_implicit_2: $(srcs-IB_implicit_2.o)
	$(CLINKER) -o $@ $(srcs-IB_implicit_2.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IBFE_explicit_0: $(srcs-IBFE_explicit_0.o)
	$(CLINKER) -o $@ $(srcs-IBFE_explicit_0.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IBFE_explicit_1: $(srcs-IBFE_explicit_1.o)
	$(CLINKER) -o $@ $(srcs-IBFE_explicit_1.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IBFE_explicit_2: $(srcs-IBFE_explicit_2.o)
	$(CLINKER) -o $@ $(srcs-IBFE_explicit_2.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IBFE_explicit_3: $(srcs-IBFE_explicit_3.o)
	$(CLINKER) -o $@ $(srcs-IBFE_explicit_3.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IBFE_explicit_4: $(srcs-IBFE_explicit_4.o)
	$(CLINKER) -o $@ $(srcs-IBFE_explicit_4.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IBFE_explicit_5: $(srcs-IBFE_explicit_5.o)
	$(CLINKER) -o $@ $(srcs-IBFE_explicit_5.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IBFE_explicit_6: $(srcs-IBFE_explicit_6.o)
	$(CLINKER) -o $@ $(srcs-IBFE_explicit_6.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IBFE_explicit_7: $(srcs-IBFE_explicit_7.o)
	$(CLINKER) -o $@ $(srcs-IBFE_explicit_7.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IBFE_explicit_8: $(srcs-IBFE_explicit_8.o)
	$(CLINKER) -o $@ $(srcs-IBFE_explicit_8.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

IMP_explicit_0: $(srcs-IMP_explicit_0.o)
	$(CLINKER) -o $@ $(srcs-IMP_explicit_0.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

adv_diff_0: $(srcs-adv_diff_0.o)
	$(CLINKER) -o $@ $(srcs-adv_diff_0.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

adv_diff_1: $(srcs-adv_diff_1.o)
	$(CLINKER) -o $@ $(srcs-adv_diff_1.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

adv_diff_2: $(srcs-adv_diff_2.o)
	$(CLINKER) -o $@ $(srcs-adv_diff_2.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

advect: $(srcs-advect.o)
	$(CLINKER) -o $@ $(srcs-advect.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

navier_stokes_0: $(srcs-navier_stokes_0.o)
	$(CLINKER) -o $@ $(srcs-navier_stokes_0.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

navier_stokes_1: $(srcs-navier_stokes_1.o)
	$(CLINKER) -o $@ $(srcs-navier_stokes_1.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

navier_stokes_2: $(srcs-navier_stokes_2.o)
	$(CLINKER) -o $@ $(srcs-navier_stokes_2.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

navier_stokes_3: $(srcs-navier_stokes_3.o)
	$(CLINKER) -o $@ $(srcs-navier_stokes_3.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

navier_stokes_4: $(srcs-navier_stokes_4.o)
	$(CLINKER) -o $@ $(srcs-navier_stokes_4.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

navier_stokes_5: $(srcs-navier_stokes_5.o)
	$(CLINKER) -o $@ $(srcs-navier_stokes_5.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

navier_stokes_6: $(srcs-navier_stokes_6.o)
	$(CLINKER) -o $@ $(srcs-navier_stokes_6.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

ConstraintIB_eel2d: $(srcs-ConstraintIB_eel2d.o)
	$(CLINKER) -o $@ $(srcs-ConstraintIB_eel2d.o) $(IBAMR_LDFLAGS) $(IBAMR_LIB) $(PETSC_LIB)

$(libibamr_libname) : $(srcs.o) | $$(@D)/.DIR
	$(call quiet,CLINKER) $(sl_linker_args) -o $@ $^ $(PETSC_EXTERNAL_LIB_BASIC)
ifneq ($(DSYMUTIL),true)
	$(call quiet,DSYMUTIL) $@
endif

$(libibamr_static) : obj := $(srcs.o)

define ARCHIVE_RECIPE_WIN32FE_LIB
  @$(RM) $@ $@.args
  @cygpath -w $^ > $@.args
  $(call quiet,AR) $(AR_FLAGS) $@ @$@.args
  @$(RM) $@.args
endef

define ARCHIVE_RECIPE_DEFAULT
  @$(RM) $@
  $(call quiet,AR) $(AR_FLAGS) $@ $^
  $(call quiet,RANLIB) $@
endef

%.$(AR_LIB_SUFFIX) : $$(obj) | $$(@D)/.DIR
	$(if $(findstring win32fe lib,$(AR)),$(ARCHIVE_RECIPE_WIN32FE_LIB),$(ARCHIVE_RECIPE_DEFAULT))

# The package libraries technically depend on each other (not just in an order-only way), but only
# ABI changes like new or removed symbols requires relinking the dependent libraries.  ABI should
# only occur when a header is changed, which would trigger recompilation and relinking anyway.
# RELINK=1 causes dependent libraries to be relinked anyway.
ifeq ($(RELINK),1)
  libdep_true = $$(libdep)
  libdep_order =
else
  libdep_true =
  libdep_order = $$(libdep)
endif
$(libpetscpkgs_libname) : $(libdep_true) | $(libdep_order) $$(@D)/.DIR
	$(call quiet,CLINKER) $(sl_linker_args) -o $@ $^ $(PETSC_EXTERNAL_LIB_BASIC)
ifneq ($(DSYMUTIL),true)
	$(call quiet,DSYMUTIL) $@
endif

%.$(SL_LINKER_SUFFIX) : $(call libname_function,%)
	@ln -sf $(notdir $<) $@

$(call soname_function,%) : $(call libname_function,%)
	@ln -sf $(notdir $<) $@

$(OBJDIR)/%.o : %.c | $$(@D)/.DIR
	$(IBAMR_COMPILE.c) $(abspath $<) -o $@

$(OBJDIR)/%.o : %.cpp | $$(@D)/.DIR
	$(IBAMR_COMPILE.cxx) $(abspath $<) -o $@

$(OBJDIR)/%.o : %.cu | $$(@D)/.DIR
	$(IBAMR_COMPILE.cu) $(abspath $<) -o $@ # Compile first so that if there is an error, it comes from a normal compile
	@$(IBAMR_GENDEPS.cu) $(abspath $<) -o $(@:%.o=%.d) # Generate the dependencies for later

# We need to do this nonsense because the Fortran preprocessor only operates in "traditional" mode, preventing concatenation
$(OBJDIR)/%.o : %.F | $$(@D)/.DIR
ifeq ($(FC_MODULE_OUTPUT_FLAG),)
	cd $(MODDIR) && $(CPP) $(FCPPFLAGS) $(abspath $<) | $(FC) -c $(FC_FLAGS) $(FFLAGS) $(FC_DEPFLAGS) -o $(abspath $@)
else
	@$(CP) $(abspath $<) $(abspath $<).C
	@$(CPP) $(FCPPFLAGS) $(abspath $<).C > $@.F
	$(IBAMR_COMPILE.F) $@.F -o $@ $(FC_MODULE_OUTPUT_FLAG)$(MODDIR)
	@$(RM) $(abspath $<).C
endif

%.F : %.f.m4
	$(M4) -DTOP_SRCDIR=$(IBAMR_DIR) -DSAMRAI_FORTDIR=$(IBAMR_DIR)/$(IBAMR_ARCH)/include $< > $@

%/.DIR :
	@mkdir -p $(@D)
	@touch $@

.PRECIOUS: %/.DIR

allobj.d := $(srcs.o:%.o=%.d)
# Tell make that allobj.d are all up to date.  Without this, the include
# below has quadratic complexity, taking more than one second for a
# do-nothing build of PETSc (much worse for larger projects)
$(allobj.d) : ;

-include $(allobj.d)
