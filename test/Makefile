# ------------------------------------------------------------ #
# Common Makefile for all scripts under src/.
#
# Authors:
#     Kee-Myoung Nam, Department of Systems Biology, Harvard Medical School
# Last updated:
#     4/11/2019
# ------------------------------------------------------------ #

CXX = g++

CXXFLAGS = --std=c++14 -frounding-math

INCLUDES = -I/home/kmn7/lib2/boost-1.68.0/build/include -I/home/kmn7/lib2/eigen-3.3.7 -I/n/groups/gunawardena/shared/CGAL-5.0.1/include 

LDFLAGS = -L/home/kmn7/lib2/boost-1.68.0/build/lib 

LIBS = -lboost_unit_test_framework -lgmp -lmpfr -lmpc

SRC_DIR = src

SRC_FILES = src/testPolynomial.cpp

BIN_DIR = bin

BIN_FILES = $(addprefix $(BIN_DIR)/, $(notdir $(patsubst %.cpp, %.out, $(SRC_FILES))))

all: $(BIN_FILES)

$(BIN_DIR)/%.out: $(SRC_DIR)/%.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $< -o $(@:.out=) $(LIBS)

clean: $(patsubst %.out, %.clean, $(BIN_FILES))

%.clean:
	$(RM) $(@:.clean=)
