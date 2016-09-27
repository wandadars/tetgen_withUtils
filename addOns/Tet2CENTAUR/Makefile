
SRC1 = Tet_2_CENTAUR.c
TARGET1 = Tet_2_CENTAUR
HOME = /home/neal
BIN = bin
CC = gcc
FLAGS = -lm

$(TARGET1) : $(SRC1)
	$(CC) $(FLAGS) $(SRC1) -o $(TARGET1)

install :
	-mv $(TARGET1) $(HOME)/$(BIN)

all : $(TARGET1) 

clean :
	-rm $(TARGET1)  *.o core
 
