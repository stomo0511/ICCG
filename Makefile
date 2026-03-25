CXX = g++
CXXFLAGS = -Wno-nan-infinity-disabled -O3 -march=native -ffast-math -std=c++17

MKL_INCLUDE = $(MKLROOT)/include
MKL_LIB = $(MKLROOT)/lib/intel64
MKL_LIBS = -lmkl_rt -lpthread -lm -ldl

UNAME = $(shell uname)
ifeq ($(UNAME), Linux)
	CXXFLAGS := -fopenmp -DUSEMKL $(CXXFLAGS) -I$(MKL_INCLUDE)
	LDFLAGS = -L$(MKL_LIB) $(MKL_LIBS)
endif

BIN_NOPRE := cg_crs
BIN_JAC   := dcg_crs
BIN_IC    := iccg_crs
BIN_PIC   := piccg_crs
BIN_ABMC  := abmc_crs

TARGET = $(BIN_NOPRE) $(BIN_JAC) $(BIN_IC) $(BIN_PIC) $(BIN_ABMC)

HDRS := crs_io.hpp precond.hpp
SRCS := crs_io.cpp cg_crs.cpp

IHDRS := $(HDRS) ic0.hpp
ISRCS := $(SRCS) ic0.cpp

PIHDRS := $(IHDRS) color.hpp
PISRCS := $(ISRCS) color.cpp

ABMCHDRS := $(PIHDRS) block.hpp
ABMCSRCS := $(PISRCS) block.cpp

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

clean:
	rm -f $(TARGET)