
INC	:= -I. -I..
LIB	:= -lm

GBM:	GBM.c Makefile
	icc -std=c99 -qopenmp -mkl -qopt-report GBM.c -o GBM $(INC) $(LIB)

clean:
	rm -f GBM

