SRC_DIR=src
OBJ_DIR=obj

cc=gcc
CFLAGS=-O3 -I(HEADER_DIR)
OPENMPFLAGS=-fopenmp

OBJ=$(OBJ_DIR)/sample.o


all: $(OBJ_DIR)

$(OBJ_DIR):
	mkdir $(OBJ_DIR)

clean:
	rm -f $(OBJ)
