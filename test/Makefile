SRC = main.f90
OBJ = $(SRC:.c=.o)

FC = gfortran

PROG = test

LIBS=foo

prog: $(OBJ) ../lib$(LIBS).a
	$(FC) $(OBJ) -L../ -l$(LIBS) -o $(PROG)

clean:
	rm -f *.o $(PROG) *.mod
