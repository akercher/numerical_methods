# -----------------------
# Andrew Kercher
# CSI 703: Assignment #6
# 03/30/2010
# -----------------------

# EXE = qtree
# OBJ = quadtree.o qmain.o 
# LIB = -lglut -lGLU -lGL
# CMP = gcc -c
# LNK = gcc

# $(EXE): $(OBJ)
# 	$(LNK) -o $@ $(OBJ) $(LIB)

# clean:
# 	rm -f *.o
# 	rm -f $(EXE)

# %.o: %.c
# 	$(CMP) $?

#DAT = circle_surf_no_uniform_test.zfem
DAT = circle_surf_no_uniform.zfem
EXE = akercher_assignment06
OBJ = quadtreeP.o akercher_assignment06.o
SRC = akercher_assignment06.c
LIB = -lglut -lGLU -lGL
CMP = gcc -c
LNK = gcc

all: $(EXE) $(DAT)
	./$?

$(EXE): $(OBJ)
	$(LNK) -o $@ $(OBJ) $(LIB)

%.o: %.c
	$(CMP) $?
clean: 
	rm -f $(EXE)
	rm -f *.o
	rm -f *~
