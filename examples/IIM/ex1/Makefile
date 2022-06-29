######################################################################
## Here specify the location of the IBAMR source and the location
## where IBAMR has been built.

IBAMR_SRC_DIR = /Users/mohammad/sfw/Capitan/ibamr/source/IBAMR-SHARP
IBAMR_BUILD_DIR = /Users/mohammad/sfw/Capitan/ibamr/aminkolah-sharp-objs-opt

######################################################################
## Include variables specific to the particular IBAMR build.
include $(IBAMR_BUILD_DIR)/config/make.inc

## Needed for Xcode to capture compiler errors and warnings.
ifdef XCODE_VERSION_ACTUAL
CXXFLAGS += -fno-color-diagnostics
endif

######################################################################
## Build the application.
##
## NOTE: The following assumes that all .cpp files in the present
##       directory are used to build the executable.
SRC = $(wildcard *.cpp)
CPPFLAGS += -MD -MP
PDIM = 2
OBJS = $(SRC:%.cpp=%.o) $(IBAMR_LIB_2D) $(IBTK_LIB_2D)
main3d_thin: $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -DNDIM=$(PDIM) -o main2d


clean:
	$(RM) main2d_thin exe
	$(RM) *.o *.lo *.objs *.ii *.int.c 
	$(RM) -r .libs

-include $(SRC:%.cpp=%.d)
