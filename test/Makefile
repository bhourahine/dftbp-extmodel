SRC = main.f90
OBJ = $(SRC:.c=.o)

PROG = test

LIBS=foo

prog: $(OBJ) ../lib$(LIBS).a
	$(FC) $(FFLAGS) $(OBJ) -L../ -l$(LIBS) -o $(PROG)

clean:
	rm -f *.o $(PROG) *.mod
