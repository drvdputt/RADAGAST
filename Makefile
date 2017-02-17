CPP=g++

INCDIR=./include

SRCDIR=./src
SOURCES=$(wildcard $(SRCDIR)/*cpp)

OBJDIR=../obj
OBJECTS=$(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SOURCES))

PROGRAM=../bin/test

CPPFLAGS=-I$(INCDIR) -O3 -std=c++11 -Wall -Wextra

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CPP) -o $(@) $^

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CPP) $(CPPFLAGS) -c -o $@ $< 

clean:
	rm $(OBJDIR)/*.o $(PROGRAM)
