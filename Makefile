SRC_DIR=src
OBJ_DIR=obj

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

test: $(OBJ) main.cpp
	$(cc) $(CFLAGS) $(OBJ) main.cpp -o test.out

call: $(OBJ) BS_call.cpp
	$(cc) $(CFLAGS) $(OBJ) BS_call.cpp -o call.out

lookback: $(OBJ) BS_lookback.cpp
	$(cc) $(CFLAGS) $(OBJ) BS_lookback.cpp -o lookback.out

barrier: $(OBJ) BS_barrier.cpp
	$(cc) $(CFLAGS) $(OBJ) BS_barrier.cpp -o barrier.out

compound: $(OBJ) compound.cpp
	$(cc) $(CFLAGS) $(OBJ) compound.cpp -o compound.out

clean:
	rm -f $(OBJ) *.out
