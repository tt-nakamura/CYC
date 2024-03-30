NTL = -lntl -lgmp
OBJ = CYC.o Period.o CYCFactoring.o ZZFactoring.o mpqs.o rho.o ZZlib.o

example: example.o $(OBJ)
	g++ example.o $(OBJ) $(NTL)
table1: table1.o $(OBJ)
	g++ table1.o $(OBJ) $(NTL)