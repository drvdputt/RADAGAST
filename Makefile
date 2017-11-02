OS:=$(shell uname)

ifeq ($(OS), Darwin)
	CXX=clang++
else
	CXX=g++
endif

# Binary target
BINDIR=../bin
$(shell mkdir -p $(BINDIR))
PROGRAM=$(BINDIR)/test

# Main source
SRCDIR=./src
SOURCES=$(wildcard $(SRCDIR)/*cpp)
HEADERS=$(wildcard $(SRCDIR)/*h)

# Includes
INCDIR=./include
EIGENDIR=./eigen3
GSLDIR=$(HOME)/.local/gsl/include

# Build directory: objects and dependency files
OBJDIR=../obj
$(shell mkdir -p $(OBJDIR) >/dev/null)
OBJECTS=$(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SOURCES))

# Archive the objects into a library
LIBDIR=../lib
$(shell mkdir -p $(LIBDIR) >/dev/null)
LIBFILE=$(LIBDIR)/libgasmodule.a
# The archiver used
AR=ar

# Automatic dependency generation
DEPENDS=$(OBJECTS:.o=.d)
$(info $$DEPENDS are [$(DEPENDS)])
# Flags that tell to compiler to create dependency files. Each dependencency file is actually
# contains instructions which fit in a makefile, and they will be included later. We create them
# here with a different name to prevent some weird conditions from occuring. They are renamed after
# generation.
DEPFLAGS=-MT $@ -MMD -MP -MF $(OBJDIR)/$*.Td

# Compile flags
CXXFLAGS=$(DEPFLAGS) -I$(INCDIR) -isystem$(GSLDIR) -isystem$(EIGENDIR) -O0 -g \
-std=c++14 -Wall -Wextra -Werror=return-type -pedantic -DREPOROOT=\""$(shell pwd)"\"

# The final target
all: $(PROGRAM) $(LIBFILE)

# Linking step
$(PROGRAM): $(OBJECTS)
	$(CXX) -o $@ $^

$(LIBFILE): $(OBJECTS)
	$(AR) rcs $@ $^

# Compiling step + dependency generation
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(OBJDIR)/%.d
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	@mv -f $(OBJDIR)/$*.Td $(OBJDIR)/$*.d

# Include the extra rules for all the objects files provided in the .d files.
# They will now depend on the correct headers.
include $(DEPENDS)

# Do nothing if the .d files don't exist yet
$(OBJDIR)/%.d: ;

# Never automatically delete .d files
.PRECIOUS: $(OBJDIR)/%.d

clean:
	rm $(OBJDIR)/*.o $(OBJDIR)/*.d $(OBJDIR)/*.Td $(PROGRAM)
