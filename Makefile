.PHONY: release
# Make should be run with the git repo as the current working directory

# Compiler
# ========

OS:=$(shell uname)

ifeq ($(OS), Darwin)
	CXX=clang++
else
	CXX=g++
endif

# The archiver used
AR=ar

# Options, to be customized by user.
# Debug build
OPTFLAGS=-O0 -g -Wall -Wextra -Wno-missing-braces -Wno-sign-compare -Werror=return-type -pedantic #-Wconversion
# Release build
release: OPTFLAGS=-O3 -Wall -Wextra -Wno-missing-braces -Werror=return-type -pedantic -DSILENT -Wno-sign-compare

# Target directories (this group of directories might be relocatable, but I have never tested
# this)
# ==================

PREFIX?=$(shell pwd)/..

# Binary target
BINDIR=$(PREFIX)/bin

# Objects and dependency files
OBJDIR=$(PREFIX)/obj

# Archive the objects that will be created into a library, for easy linking by the client code.
LIBDIR=$(PREFIX)/lib

# Symlink the headers describing the public interface, to have a clear and minimal include path
# for the client code.
PUBLIC_INCLUDEDIR=$(PREFIX)/include

# Source directories
# ==================

# Core source code
SRCDIR=./src
COREDIR=$(SRCDIR)/core
MAINDIR=$(SRCDIR)/mains

# Includes
INCDIR=./include
EIGENDIR=./eigen3
GSLDIR=/usr/lib/x86_64-linux-gnu

# Include flags. Adjust the different paths above.
INCFLAGS=-I$(COREDIR) -I$(INCDIR) -I$(GSLDIR)/include -isystem$(EIGENDIR)

# Linker flags (not that leaving out -lgslcblas will link, but it's a lazy link an will crash at
# runtime)
LNKFLAGS=-L$(GSLDIR)/lib -lgsl -lgslcblas -lm

# Source files
# ============

# Core source
SOURCES=$(wildcard $(COREDIR)/*cpp)
HEADERS=$(wildcard $(COREDIR)/*h)

# Executable main sources
MAINS=$(wildcard $(MAINDIR)/*cpp)

# Target files
# ============

# All the targets will be placed in the parent directory of the repo.

# Make the necessary subdirectories
OBJ_SUBDIRS=$(patsubst $(SRCDIR)/%, $(OBJDIR)/%, $(MAINDIR) $(COREDIR))
$(shell mkdir -p $(BINDIR) $(OBJ_SUBDIRS) $(LIBDIR))

# Binary
BINARIES=$(patsubst $(MAINDIR)/%.cpp, $(BINDIR)/%, $(MAINS))

# Objects
COREOBJECTS=$(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SOURCES))
MAINOBJECTS=$(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(MAINS))

# Library
LIBNAME=gasmodule
CORELIB=$(LIBDIR)/lib$(LIBNAME).a

# We also make an include directory with symlinks to some headers in the repo. These should be
# the only headers that a client will ever need to include.
$(shell mkdir -p $(PUBLIC_INCLUDEDIR) >/dev/null)
INTERFACEHEADERS=GasInterface.h GasState.h GrainInterface.h
INTERFACELINKS=$(patsubst %, $(PUBLIC_INCLUDEDIR)/%, $(INTERFACEHEADERS))

# Dependency files.
DEPENDS=$(COREOBJECTS:.o=.d)

# Flags that tell to compiler to create dependency files. Each dependencency file is actually
# contains instructions which fit in a makefile, and they will be included later. We create them
# here with a different name to prevent some weird conditions from occuring. They are renamed
# after generation.
DEPFLAGS=-MT $@ -MMD -MP -MF $(OBJDIR)/$*.Td

# Build commands and dependencies
# ===============================
# Cheatsheet:
# $^ == all prerequisites
# $< == leftmost prerequisite
# $@ == target

# All the flags. Do not change. Change the variables above instead.
CXXFLAGS=$(OPTFLAGS) $(INCFLAGS) $(DEPFLAGS) -std=c++14 -DREPOROOT=\""$(shell pwd)"\"

# The final targets
all release: $(CORELIB) $(INTERFACELINKS) $(BINARIES)

# Linking step
# Make an exception for tests, as these need the list of objects
$(BINDIR)/test: $(BINDIR)/%: $(OBJDIR)/mains/%.o $(COREOBJECTS)
	$(CXX) -o $@ $^ $(LNKFLAGS)
# The rest of the binaries can be linked to the library file
$(BINDIR)/%: $(OBJDIR)/mains/%.o $(CORELIB)
	$(CXX) -o $@ $< -L$(LIBDIR) -l$(LIBNAME) $(LNKFLAGS)

# Create archive
$(CORELIB): $(COREOBJECTS)
	$(AR) rcs $@ $^

# Create symlinks to public interface headers (ln -s [target] [link name])
$(INTERFACELINKS): $(PUBLIC_INCLUDEDIR)/%: $(COREDIR)/%
	ln -sf "$(shell pwd)"/$< $@

# Compiling step + dependency generation
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(OBJDIR)/%.d
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
	@mv -f $(OBJDIR)/$*.Td $(OBJDIR)/$*.d

# Include the extra rules for all the objects files provided in the .d files. They will now
# depend on the correct headers.
include $(DEPENDS)

# Do nothing if the .d files don't exist yet
$(OBJDIR)/%.d: ;

# Never automatically delete .d files
.PRECIOUS: $(OBJDIR)/%.d

clean:
	rm -rf $(OBJDIR)/* $(BINDIR)/* $(PUBLIC_INCLUDEDIR)/* $(LIBDIR)/*

.PHONY: doc
doc:
	(cd dox && doxygen doxygen.conf)
