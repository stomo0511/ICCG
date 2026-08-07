CXX := g++
CXXFLAGS := -O3 -march=native -std=c++17
LDFLAGS := 

UNAME = $(shell uname)
USER_NAME = $(shell whoami)
HOST_NAME = $(shell hostname)

ifeq ($(UNAME), Darwin)  # MacOS
	CXXFLAGS := -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include $(CXXFLAGS)
	LDFLAGS := -L/opt/homebrew/opt/libomp/lib -lomp
else ifeq ($(UNAME), Linux)
	MKL_INCLUDE = $(MKLROOT)/include
	MKL_LIB = $(MKLROOT)/lib/intel64
	MKL_LIBS = -lmkl_rt -lpthread -lm -ldl

	CXXFLAGS := -fopenmp -DUSEMKL $(CXXFLAGS) -I$(MKL_INCLUDE)
	LDFLAGS := -L$(MKL_LIB) $(MKL_LIBS)

	ifeq ($(USER_NAME), a30029)  # 北大
		CXX := icpc
	endif
endif

BIN_NOPRE := cg
BIN_JAC   := dcg
BIN_IC    := iccg
BIN_PIC   := mciccg
BIN_ABMC  := bmciccg
BIN_LLT   := llt_eval
BIN_LLT_IC    := llt_ic
BIN_LLT_MCIC  := llt_mcic
BIN_LLT_BMCIC := llt_bmcic

TARGET = $(BIN_NOPRE) $(BIN_JAC) $(BIN_IC) $(BIN_PIC) $(BIN_ABMC)

HDRS := crs_io.hpp precond.hpp
SRCS := crs_io.cpp cg_crs.cpp

IHDRS := $(HDRS) ic0.hpp
ISRCS := $(SRCS) ic0.cpp

PIHDRS := $(IHDRS) color.hpp
PISRCS := $(ISRCS) color.cpp

ABMCHDRS := $(PIHDRS) block.hpp
ABMCSRCS := $(PISRCS) block.cpp

LLTHDRS := crs_io.hpp color.hpp block.hpp ic0.hpp
LLTSRCS := crs_io.cpp color.cpp block.cpp ic0.cpp llt_eval.cpp
LLT_ERR_HDRS := crs_io.hpp ic0.hpp llt_err.hpp
LLT_IC_SRCS := crs_io.cpp ic0.cpp llt_err.cpp llt_ic.cpp
LLT_MCIC_SRCS := crs_io.cpp color.cpp ic0.cpp llt_err.cpp llt_mcic.cpp
LLT_BMCIC_SRCS := crs_io.cpp color.cpp block.cpp ic0.cpp llt_err.cpp llt_bmcic.cpp
LLT_CXXFLAGS := $(filter-out -DUSEMKL,$(CXXFLAGS))
LLT_LDFLAGS := $(LDFLAGS)
ifeq ($(UNAME), Linux)
	LLT_LDFLAGS :=
endif

all: $(TARGET)

llt: $(BIN_LLT_IC) $(BIN_LLT_MCIC) $(BIN_LLT_BMCIC)

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

$(BIN_LLT): $(LLTHDRS) $(LLTSRCS)
	$(CXX) $(LLT_CXXFLAGS) -DIC -DABMC -o $@ $(LLTSRCS) $(LLT_LDFLAGS)

$(BIN_LLT_IC): $(LLT_ERR_HDRS) $(LLT_IC_SRCS)
	$(CXX) $(LLT_CXXFLAGS) -DIC -o $@ $(LLT_IC_SRCS) $(LLT_LDFLAGS)

$(BIN_LLT_MCIC): $(LLT_ERR_HDRS) color.hpp $(LLT_MCIC_SRCS)
	$(CXX) $(LLT_CXXFLAGS) -DIC -o $@ $(LLT_MCIC_SRCS) $(LLT_LDFLAGS)

$(BIN_LLT_BMCIC): $(LLT_ERR_HDRS) color.hpp block.hpp $(LLT_BMCIC_SRCS)
	$(CXX) $(LLT_CXXFLAGS) -DIC -DABMC -o $@ $(LLT_BMCIC_SRCS) $(LLT_LDFLAGS)

clean:
	rm -f $(TARGET) $(BIN_LLT) $(BIN_LLT_IC) $(BIN_LLT_MCIC) $(BIN_LLT_BMCIC)
