SRC_DIR=src
OBJ_DIR=obj
EXE_DIR=executable
SIM_DIR=simulations

cc=g++
CFLAGS=-O3 -I$(SRC_DIR)
OPENMPFLAGS=-fopenmp

SRC= structural_params.cpp \
	 models.cpp \
	 estimator.cpp

OBJ=$(OBJ_DIR)/structural_params.o \
	$(OBJ_DIR)/models.o \
	$(OBJ_DIR)/estimator.o

all: $(OBJ_DIR) $(EXE_DIR) call lookback barrier compound

$(OBJ_DIR):
	mkdir $(OBJ_DIR)

$(EXE_DIR):
	mkdir $(EXE_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(cc) $(CFLAGS) -c $< -o $@

call: $(OBJ) $(SIM_DIR)/BS_call.cpp
	$(cc) $(CFLAGS) $(OBJ) $(SIM_DIR)/BS_call.cpp -o $(EXE_DIR)/call.out

lookback: $(OBJ) $(SIM_DIR)/BS_lookback.cpp 
	$(cc) $(CFLAGS) $(OBJ) $(SIM_DIR)/BS_lookback.cpp -o $(EXE_DIR)/lookback.out

barrier: $(OBJ) $(SIM_DIR)/BS_barrier.cpp
	$(cc) $(CFLAGS) $(OBJ) $(SIM_DIR)/BS_barrier.cpp -o $(EXE_DIR)/barrier.out

compound: $(OBJ) $(SIM_DIR)/compound.cpp
	$(cc) $(CFLAGS) $(OBJ) $(SIM_DIR)/compound.cpp -o $(EXE_DIR)/compound.out

clean:
	rm -f $(OBJ) *.out *.txt executable/*
	rmdir $(OBJ_DIR)
	rmdir $(EXE_DIR)
