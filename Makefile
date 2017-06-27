CPP=g++

# Binary target
PROGRAM=../bin/test

# Main source
SRCDIR=./src
SOURCES=$(wildcard $(SRCDIR)/*cpp)
HEADERS=$(wildcard $(SRCDIR)/*h)

# Includes
INCDIR=./include
EIGENDIR=./eigen3

# Build directory: objects and dependency files
OBJDIR=../obj
$(shell mkdir -p $(OBJDIR) >/dev/null)
OBJECTS=$(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SOURCES))

# Automatic dependency generation
DEPENDS=$(OBJECTS:.o=.d)
$(info $$DEPENDS are [$(DEPENDS)])
DEPFLAGS=-MT $@ -MMD -MP -MF $(OBJDIR)/$*.Td

# Compile flags
CPPFLAGS=$(DEPFLAGS) -I$(INCDIR) -isystem$(EIGENDIR) -O0 -g -std=c++14 -Wall\
-Wextra -Werror=return-type -pedantic -DREPOROOT=\""$(shell pwd)"\"

# The final target
all: $(PROGRAM)

# Linking step
$(PROGRAM): $(OBJECTS)
	$(CPP) -o $(@) $^

# Compiling step + dependency generation
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(OBJDIR)/%.d
	$(CPP) $(CPPFLAGS) -c -o $@ $<
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
