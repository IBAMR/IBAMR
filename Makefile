include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables
include $(IBAMR_DIR)/$(IBAMR_ARCH)/lib/petsc/conf/ibamrvariables

IBAMR_ARCH := $(if $(IBAMR_ARCH),$(IBAMR_ARCH),build)

all : $(IBAMR_ARCH)/conf/configure.log $(IBAMR_ARCH)/gmakefile
	$(MAKE) -C $(IBAMR_ARCH) -f gmakefile
	@echo "Build complete in $(IBAMR_ARCH).  Use make test to test."

$(IBAMR_ARCH)/conf/configure.log:
	./configure.new --IBAMR_ARCH=$(IBAMR_ARCH) --download-muparser --download-silo --download-samrai --with-mpi-dir=$(PETSC_DIR)/$(PETSC_ARCH) --with-hdf5-dir=$(PETSC_DIR)/$(PETSC_ARCH) --with-eigen-dir=$(PETSC_DIR)/$(PETSC_ARCH)

$(IBAMR_ARCH)/gmakefile: ./config/gmakegen.py
	$(MAKE) -C $(IBAMR_ARCH) -f bootstrap.mk gmakefile

test : all
	$(MAKE) -C $(IBAMR_ARCH) test

clean :
	$(MAKE) -C $(IBAMR_ARCH) clean

.PHONY: all test clean
