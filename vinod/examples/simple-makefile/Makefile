## ---------------------------------------------------------------------
##
## Copyright (c) 2025 - 2025 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

# Simple Makefile for standalone IBAMR example
# This assumes IBAMR is installed and available via pkg-config or environment variables

# Executable name
TARGET = simple-makefile-example

# Source files
SOURCES = main.cpp

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Compiler
CXX ?= mpic++

# Try to detect IBAMR installation
IBAMR_PREFIX ?= /usr/local

# Include directories
INCLUDES = -I$(IBAMR_PREFIX)/include \
           -I$(IBAMR_PREFIX)/include/ibamr \
           -I$(IBAMR_PREFIX)/include/ibtk

# Library directories and libraries
LDFLAGS = -L$(IBAMR_PREFIX)/lib
LIBS = -libamr2d -libtk2d -lsamrai

# Compiler flags
CXXFLAGS = -std=c++17 -O2 -Wall -DNDIM=2 $(INCLUDES)

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

# Compile source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJECTS) $(TARGET)

# Run the example
run: $(TARGET)
	mpirun -np 1 ./$(TARGET)

# Run with multiple processors
run-parallel: $(TARGET)
	mpirun -np 4 ./$(TARGET)

.PHONY: all clean run run-parallel
