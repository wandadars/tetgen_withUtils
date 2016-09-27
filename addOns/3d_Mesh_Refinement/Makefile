
SRC1 = 3d_amr_code.c
TARGET1 = 3d_amr_code
HOME = /home/neal
BIN = bin
CC = gcc
FLAGS = -lm

$(TARGET1) : $(SRC1)
	-$(CC) $(FLAGS) $(SRC1) -o $(TARGET1)

install :
	-mv $(TARGET1) $(HOME)/$(BIN)

all : $(TARGET1) 

clean :
	-rm $(TARGET1)  *.o core
 
