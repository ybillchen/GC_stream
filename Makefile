# BSD 3-Clause License
# Copyright (c) 2023 Yingtian Chen
# All rights reserved.

LINK = $(CXX)

AGAMA_PATH = /Users/ybchen/opt/anaconda3/lib/python3.8/site-packages/agama
CXXFLAGS   += -I$(AGAMA_PATH)/src -fPIC -fopenmp -Wall -O2 -march=native -std=c++11

# needed since agama relies on python lib
EXE_FLAGS  = -Wl,-rpath,/Users/ybchen/opt/anaconda3/lib /Users/ybchen/opt/anaconda3/lib/libpython3.8.dylib

# set up folder names
SRCDIR   = GC_stream
OBJDIR   = obj
EXEDIR   = exe
TESTSDIR = tests

# sources of the main library
SOURCES  = potential_varying.cpp \


# test and example programs
TESTSRCS = test_demo.cpp \
		   test_orbit.cpp \

LIBNAME  = agamaVarying.so
OBJECTS  = $(patsubst %.cpp,$(OBJDIR)/%.o,$(SOURCES))
TESTEXE  = $(patsubst %.cpp,$(EXEDIR)/%.exe,$(TESTSRCS))
CXXFLAGS += -I$(SRCDIR)

# this is the default target (build all), if make is launched without parameters
all:  $(LIBNAME) $(TESTEXE)

# one may recompile just the shared library by running 'make lib'
lib:  $(LIBNAME)

# one may recompile just the test files by running 'make test'
test:  $(TESTEXE)

# the .so format is hard to deal with in Mac. So let's just create a symlink here
$(LIBNAME):  $(OBJECTS) Makefile
	@[ -f agama.so -a ! -L agama.so ] && rm agama.so || true
	@[ -L agama.so ] || ln -s $(AGAMA_PATH)/agama.so agama.so
	$(LINK) -shared -o $(LIBNAME) $(OBJECTS) agama.so $(CXXFLAGS)

# the same: create symlink for .so files
$(EXEDIR)/%.exe:  $(TESTSDIR)/%.cpp $(LIBNAME)
	@mkdir -p $(EXEDIR)
	@[ -f $(EXEDIR)/$(LIBNAME) -a ! -L $(EXEDIR)/$(LIBNAME) ] && rm $(EXEDIR)/$(LIBNAME) || true
	@[ -L $(EXEDIR)/$(LIBNAME) ] || ln -s ../$(LIBNAME) $(EXEDIR)/$(LIBNAME)
	@[ -f $(EXEDIR)/agama.so -a ! -L $(EXEDIR)/agama.so ] && rm $(EXEDIR)/agama.so || true
	@[ -L $(EXEDIR)/agama.so ] || ln -s $(AGAMA_PATH)/agama.so $(EXEDIR)/agama.so
	$(LINK) -o $@ $< $(LIBNAME) agama.so $(CXXFLAGS) $(EXE_FLAGS)

$(OBJDIR)/%.o:  $(SRCDIR)/%.cpp
	@mkdir -p $(OBJDIR)
	$(CXX) -c $(CXXFLAGS) $(COMPILE_FLAGS) -o $@ $<

clean:
	rm -f *.so $(OBJDIR)/* $(EXEDIR)/*

.PHONY: clean test lib