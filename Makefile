CPP=g++

INCDIR=./include
EIGENDIR=./eigen3

SRCDIR=./src
SOURCES=$(wildcard $(SRCDIR)/*cpp)
HEADERS=$(wildcard $(SRCDIR)/*h)

SOURCESDOTH=$(patsubst %.cpp, %.h, $(SOURCES))
STANDALONEHEADERS=$(filter-out $(SOURCESDOTH), $(HEADERS))

OBJDIR=../obj
OBJECTS=$(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SOURCES))

PROGRAM=../bin/test

CPPFLAGS=-I$(INCDIR) -I$(EIGENDIR) -O0 -g -std=c++14 -Wall -Wextra -Werror=return-type

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CPP) -o $(@) $^

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(STANDALONEHEADERS)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

clean:
	rm $(OBJDIR)/*.o $(PROGRAM)
