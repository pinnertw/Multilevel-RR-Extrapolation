SRC_DIR=src
OBJ_DIR=obj

cc=g++
CFLAGS=-O3 -I$(SRC_DIR)
OPENMPFLAGS=-fopenmp

SRC= structural_params.cpp \
	 models.cpp

OBJ=$(OBJ_DIR)/structural_params.o \
	$(OBJ_DIR)/models.o

all: $(OBJ_DIR) test

$(OBJ_DIR):
	mkdir $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(cc) $(CFLAGS) -c $< -o $@

test: $(OBJ) main.cpp
	$(cc) $(CFLAGS) $(OBJ) main.cpp -o test.out

clean:
	rm -f $(OBJ) test.out
