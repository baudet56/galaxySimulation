# Author : author
# Copyright ENSICAEN
# Date: date

CC=gcc
CCFLAGS=-Wall -Wextra -pedantic -I./include -g -O3 -lm -pthread 

BIN=./bin
OBJ=./obj
SRC=./src

SRC_FILE = $(wildcard $(SRC)/*.c)
OBJ_FILE = $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SRC_FILE))

DIRECTORY=bin obj include src
EXEC=main

all: prepare $(EXEC)

prepare: | $(DIRECTORY)

$(DIRECTORY): 
	mkdir -p $@
	touch README.md

$(EXEC): $(OBJ_FILE)
	$(CC) -o $(BIN)/$@ $^ $(CCFLAGS)

$(OBJ)/imageGenerator.o: $(SRC)/imageGenerator.c ./include/imageGenerator.h ./include/star.h
	$(CC) -o $@ -c $< $(CCFLAGS)

$(OBJ)/saveData.o: $(SRC)/saveData.c ./include/saveData.h ./include/star.h
	$(CC) -o $@ -c $< $(CCFLAGS)

$(OBJ)/star.o: $(SRC)/star.c ./include/star.h
	$(CC) -o $@ -c $< $(CCFLAGS)


$(OBJ)/main.o: $(SRC)/main.c ./include/imageGenerator.h ./include/star.h ./include/saveData.h
	$(CC) -o $@ -c $< $(CCFLAGS)



clean:
	rm -rf $(OBJ) $(LIB)

properclean: clean
	rm -rf $(BIN)

archive: properclean
	tar -cvf "$(EXEC)_author.tar.gz" *
	find . -type f > file.txt
	find ./include -type f >> file.txt
	find ./src -type f >> file.txt

.PHONY: all prepare directory file clean properclean
