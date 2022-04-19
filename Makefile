SRC_DIR=src
OBJ_DIR=obj
INCLUDE=include

cc=g++
CFLAGS=-O3 -I$(INCLUDE)
OPENMPFLAGS=-fopenmp

all: $(OBJ_DIR) params test

$(OBJ_DIR):
	mkdir $(OBJ_DIR)

params: $(SRC_DIR)/params.cpp
	$(cc) $(CFLAGS) -DPRINT $< -c -o obj/params.o

generator : $(SRC_DIR)/generator.cpp
	$(cc) $(CFLAGS) $< $(OBJ_DIR)/params.o -o test.out

test : $(SRC_DIR)/generator.cpp
	$(cc) $(CFLAGS) -DPRINT $< $(OBJ_DIR)/params.o -o test.out

clean:
	rm -f $(OBJ)
