CXX = g++
# CXXFLAGS = -Wno-nan-infinity-disabled -O3 -march=native -ffast-math -std=c++17
CXXFLAGS = -O3 -march=native -std=c++17

MKL_INCLUDE = $(MKLROOT)/include
MKL_LIB = $(MKLROOT)/lib/intel64
MKL_LIBS = -lmkl_rt -lpthread -lm -ldl

UNAME = $(shell uname)
ifeq ($(UNAME), Linux)
	CXXFLAGS := -fopenmp -DUSEMKL $(CXXFLAGS) -I$(MKL_INCLUDE)
	LDFLAGS := -L$(MKL_LIB) $(MKL_LIBS)
else
	CXXFLAGS := -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
	LDFLAGS := -L/opt/homebrew/opt/libomp/lib -lomp
endif

BIN_NOPRE := cg
BIN_JAC   := dcg
BIN_IC    := iccg
BIN_PIC   := miccg
BIN_ABMC  := bmciccg
BIN_CONDEST := condest_cg

TARGET = $(BIN_NOPRE) $(BIN_JAC) $(BIN_IC) $(BIN_PIC) $(BIN_ABMC) $(BIN_CONDEST)

HDRS := crs_io.hpp precond.hpp
SRCS := crs_io.cpp cg_crs.cpp

IHDRS := $(HDRS) ic0.hpp
ISRCS := $(SRCS) ic0.cpp

PIHDRS := $(IHDRS) color.hpp
PISRCS := $(ISRCS) color.cpp

ABMCHDRS := $(PIHDRS) block.hpp
ABMCSRCS := $(PISRCS) block.cpp

CONDEST_HDRS := crs_io.hpp
CONDEST_SRCS := crs_io.cpp condest_cg.cpp

all: $(TARGET)

$(BIN_NOPRE): $(HDRS) $(SRCS)
	$(CXX) $(CXXFLAGS) -DNOPRE -o $@ $(SRCS) $(LDFLAGS)

$(BIN_JAC): $(HDRS) $(SRCS)
	$(CXX) $(CXXFLAGS) -DJAC -o $@ $(SRCS) $(LDFLAGS)

$(BIN_IC): $(IHDRS) $(ISRCS)
	$(CXX) $(CXXFLAGS) -DIC -o $@ $(ISRCS) $(LDFLAGS)

$(BIN_PIC): $(PIHDRS) $(PISRCS)
	$(CXX) $(CXXFLAGS) -DIC -DPIC -o $@ $(PISRCS) $(LDFLAGS)

$(BIN_ABMC): $(ABMCHDRS) $(ABMCSRCS)
	$(CXX) $(CXXFLAGS) -DIC -DPIC -DABMC -o $@ $(ABMCSRCS) $(LDFLAGS)

$(BIN_CONDEST): $(CONDEST_HDRS) $(CONDEST_SRCS)
	$(CXX) $(CXXFLAGS) -o $@ $(CONDEST_SRCS) $(LDFLAGS)

clean:
	rm -f $(TARGET)
