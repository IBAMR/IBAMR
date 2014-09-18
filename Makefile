IBAMR_ARCH := $(if $(IBAMR_ARCH),$(IBAMR_ARCH),build)

all : $(IBAMR_ARCH)/conf/configure.log
	$(MAKE) -C $(IBAMR_ARCH) -f gmakefile
	@echo "Build complete in $(IBAMR_ARCH).  Use make test to test."

$(IBAMR_ARCH)/conf/configure.log:
	./configure.new --IBAMR_ARCH=$(IBAMR_ARCH)

test : all
	$(MAKE) -C $(IBAMR_ARCH) test

clean :
	$(MAKE) -C $(IBAMR_ARCH) clean

.PHONY: all test clean
