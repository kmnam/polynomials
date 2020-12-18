# ------------------------------------------------------------ #
# Common Makefile for all scripts under src/.
#
# Authors:
#     Kee-Myoung Nam, Department of Systems Biology, Harvard Medical School
# Last updated:
#     4/11/2019
# ------------------------------------------------------------ #

CXX = g++-10

CXXFLAGS = --std=c++11 -frounding-math

INCLUDES = -I/usr/local/Cellar/boost/1.71.0/include

LDFLAGS = -L/usr/local/Cellar/boost/1.71.0/lib 

LIBS = -lmpfr -lmpc -lboost_unit_test_framework

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
