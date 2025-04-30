CC=g++

OPTS= --std=c++14
INCLUDES=./external_libraries/include

BUILD_DIR=./build

.PHONY = all external


all: 
	$(MAKE) external; 
	$(CC) -I $(INCLUDES) -g --std=c++14 $(BUILD_DIR)/cargs.o studnie.cpp main.cpp -o studnie

external:
	$(CC) -I $(INCLUDES) -g -c ./external_libraries/src/cargs.c -o $(BUILD_DIR)/cargs.o