# -----------------------
# Andrew Kercher
# CSI 703: Assignment #6
# 03/30/2010
# -----------------------

EXE = run
OBJ = quadtreePR.o searchCompare.o akercher_assignment06.o
SRC = akercher_assignment06.c
LIB = -lglut -lGLU -lGL
CMP = gcc -c
LNK = gcc

$(EXE): $(OBJ)
	$(LNK) -o $@ $(OBJ) $(LIB)

%.o: %.c
	$(CMP) $?
clean: 
	rm -f *.o
	rm -f $(EXE)