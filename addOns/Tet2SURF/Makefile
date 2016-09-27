
#Pass PARFLAG=1 when making to compile the parallel version. I.e. make PARFLAG=1

ifdef PARFLAG
SRC1 = Tet_2_SURF_parallel.c
TARGET1 = Tet_2_SURF
HOME = /home/neal
BIN = bin
CC = gcc
FLAGS = -lm -g -fopenmp -O3
$(info Compiling parallel Version) 
else
SRC1 = Tet_2_SURF_serial.c
TARGET1 = Tet_2_SURF
HOME = /home/neal
BIN = bin
CC = gcc
FLAGS = -lm -g -O3
$(info Compiling serial version)
endif


$(TARGET1) : $(SRC1)
	$(CC) $(FLAGS) $(SRC1) -o $(TARGET1)

install :
	-mv $(TARGET1) $(HOME)/$(BIN)

all : $(TARGET1) 

clean :
	-rm $(TARGET1)  *.o core
 
