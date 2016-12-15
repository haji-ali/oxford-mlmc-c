TARGETS = mlmc
DEBUG = -O3

all: ${TARGETS}

mlmc: mlmc.c main.c
	gcc ${DEBUG} -fPIC -Wall -o $@ $^ -lm

clean:
	${RM} *.so *.a *.o
