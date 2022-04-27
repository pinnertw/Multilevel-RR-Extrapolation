SRC_DIR=src
OBJ_DIR=obj
EXE_DIR=executable
SIM_DIR=simulations
RES_DIR=results

cc=g++
CFLAGS=-O3 -I$(SRC_DIR)
MLMC=-DMLMC

OBJ=$(OBJ_DIR)/structural_params.o $(OBJ_DIR)/models.o $(OBJ_DIR)/estimator.o

all: $(OBJ_DIR) $(EXE_DIR) $(RES_DIR) call_MLRR lookback_MLRR barrier_MLRR compound_MLRR call_MLMC lookback_MLMC barrier_MLMC compound_MLMC

$(OBJ_DIR):
	mkdir $(OBJ_DIR)

$(EXE_DIR):
	mkdir $(EXE_DIR)

$(RES_DIR):
	mkdir $(RES_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(cc) $(CFLAGS) -c $< -o $@

call_MLRR: $(OBJ) $(SIM_DIR)/BS_call.cpp
	$(cc) $(CFLAGS) $^ -o $(EXE_DIR)/$@.out

lookback_MLRR: $(OBJ) $(SIM_DIR)/BS_lookback.cpp 
	$(cc) $(CFLAGS) $^ -o $(EXE_DIR)/$@.out

barrier_MLRR: $(OBJ) $(SIM_DIR)/BS_barrier.cpp
	$(cc) $(CFLAGS) $^ -o $(EXE_DIR)/$@.out

compound_MLRR: $(OBJ) $(SIM_DIR)/compound.cpp
	$(cc) $(CFLAGS) $^ -o $(EXE_DIR)/$@.out

call_MLMC: $(OBJ) $(SIM_DIR)/BS_call.cpp
	$(cc) $(CFLAGS) $(MLMC) $^ -o $(EXE_DIR)/$@.out

lookback_MLMC: $(OBJ) $(SIM_DIR)/BS_lookback.cpp 
	$(cc) $(CFLAGS) $(MLMC) $^ -o $(EXE_DIR)/$@.out

barrier_MLMC: $(OBJ) $(SIM_DIR)/BS_barrier.cpp
	$(cc) $(CFLAGS) $(MLMC) $^ -o $(EXE_DIR)/$@.out

compound_MLMC: $(OBJ) $(SIM_DIR)/compound.cpp
	$(cc) $(CFLAGS) $(MLMC) $^ -o $(EXE_DIR)/$@.out

clean:
	rm -f $(OBJ) *.out *.txt executable/*
	rmdir $(OBJ_DIR)
	rmdir $(EXE_DIR)
