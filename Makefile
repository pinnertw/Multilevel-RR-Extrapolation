SRC_DIR=src
OBJ_DIR=obj

cc=g++
CFLAGS=-O3
OPENMPFLAGS=-fopenmp

SRC=$(SRC_DIR)/generator.cpp
OBJ=$(OBJ_DIR)/generator.o


all: $(OBJ_DIR) test #generator

$(OBJ_DIR):
	mkdir $(OBJ_DIR)

generator : $(SRC_DIR)/generator.cpp
	$(cc) $(CFLAGS) $< -o test.out

test : $(SRC_DIR)/generator.cpp
	$(cc) $(CFLAGS) $< -DPRINT -o test.out

clean:
	rm -f $(OBJ)
