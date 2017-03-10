CPP=g++

INCDIR=./include
EIGENDIR=./eigen3

SRCDIR=./src
SOURCES=$(wildcard $(SRCDIR)/*cpp)

OBJDIR=../obj
OBJECTS=$(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SOURCES))

PROGRAM=../bin/test

CPPFLAGS=-I$(INCDIR) -I$(EIGENDIR) -O3 -std=c++11 -Wall -Wextra -Werror=return-type

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CPP) -o $(@) $^

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CPP) $(CPPFLAGS) -c -o $@ $< 

clean:
	rm $(OBJDIR)/*.o $(PROGRAM)
