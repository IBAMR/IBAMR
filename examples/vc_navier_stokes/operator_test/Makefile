######################################################################
## Here specify the location of the IBAMR source and the location
## where IBAMR has been built.
IBAMR_SRC_DIR = $(HOME)/ibamr_code/IBAMR
IBAMR_BUILD_DIR = $(HOME)/ibamr_code/ibamr-objs-debug

top_srcdir = $(IBAMR_SRC_DIR)
top_builddir = $(IBAMR_BUILD_DIR)

######################################################################
## Include variables specific to the particular IBAMR build.
include $(IBAMR_BUILD_DIR)/config/make.inc

######################################################################
## Build the tester application.
SOURCES = main.C
OBJS = main.o
PDIM = 2

default: main2d

main2d: $(OBJS) $(IBAMR_LIB_2D) $(IBTK_LIB_2D)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(IBAMR_LIB_2D) $(IBTK_LIB_2D) $(LIBS) -DNDIM=$(PDIM) -o main2d

main3d: $(OBJS) $(IBAMR_LIB_3D) $(IBTK_LIB_3D)
        $(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(LIBS) -DNDIM=$(PDIM) -o main3d

clean:
	$(RM) main2d main3d
	$(RM) *.o *.lo *.objs *.ii *.int.c
	$(RM) -r .libs

