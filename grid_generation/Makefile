#------------------------------------------
# Student: Andrew Kercher
# CSI 722: Project 2, grid generation with advancing front.
# Last Updated: 12/01/2011
# Program: Makefile 
#------------------------------------------

#-----------------------
# VARIABLES: FORTRAN 90
#-----------------------

#---------------------------------
# External Modules for Fortran 90
#---------------------------------
# OBJ =  csi721_typeDef.o csi721_dataIO.o csi721_dataStruc.o csi721_sub.o csi721_project2.o 
MODF1 = csi722_typeDef.o
MODF2 = csi722_advFront.o
SRCM1 = csi722_typeDef.f90
SRCM2 = csi722_advFront.f90

#---------------------------------------------------
# Program and Output Files for csi722_project1.f90
#---------------------------------------------------
EXEF = csi722_project2.run
OBJF = csi722_project2.o 
SRCF = csi722_project2.f90

#--------------
# VARIABLES: C
#--------------
EXEC = csi722_grid.run
OBJC = csi722_gridVisual.o
SRCC = csi722_gridVisual.c

#----------------------------------------
# LIBRARIES AND COMPILER TO USE IN LINUX
#----------------------------------------
LIBF = -lblas -llapack
CMPF = gfortran -c -g 
LNKF = gfortran
LIBC = -lglut -lGLU -lGL
CMPC = gcc -c
LNKC = gcc

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

#--------------
# TARGETS: all
#--------------
all: $(EXEF)
	./$(EXEF)

#-----------------
# TARGETS: visual
#-----------------
visual: $(EXEC)
	./$(EXEC) data/faces data/pts

#---------------------------------------------
# DEPENDENCIES: Compile needed *.f90 programs
#---------------------------------------------
$(MODF1): $(SRCM1)
	$(CMPF) $?

$(MODF2): $(SRCM2)
	$(CMPF) $?

$(OBJF): $(SRCF) 
	$(CMPF) $?

$(EXEF): $(MODF1) $(MODF2) $(OBJF)
	$(LNKF) -o $@ $? $(LIBF)

#------------------------------------------------------------
# SECOND TARGET: Produce visualization using openGL and glut
#------------------------------------------------------------
$(OBJC): $(SRCC)
	$(CMPC) $?

$(EXEC): $(OBJC)
	$(LNKC) -o $@ $(OBJC) $(LIBC)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

#-------------------------------
# CLEAN EXECUTABLES AND OBJECTS
#-------------------------------
clean:
	rm -f $(EXEC)
	rm -f $(EXEF)
	rm -f *.o
	rm -f *.mod
	rm -f *~