SRC_DIR=src
OBJ_DIR=obj
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

all: $(OBJ_DIR) test

$(OBJ_DIR):
	mkdir $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(cc) $(CFLAGS) -c $< -o $@

call: $(OBJ) $(SIM_DIR)/BS_call.cpp
	$(cc) $(CFLAGS) $(OBJ) $(SIM_DIR)/BS_call.cpp -o call.out

lookback: $(OBJ) $(SIM_DIR)/BS_lookback.cpp
	$(cc) $(CFLAGS) $(OBJ) $(SIM_DIR)/BS_lookback.cpp -o lookback.out

barrier: $(OBJ) $(SIM_DIR)/BS_barrier.cpp
	$(cc) $(CFLAGS) $(OBJ) $(SIM_DIR)/BS_barrier.cpp -o barrier.out

compound: $(OBJ) compound.cpp
	$(cc) $(CFLAGS) $(OBJ) compound.cpp -o compound.out

clean:
	rm -f $(OBJ) *.out
