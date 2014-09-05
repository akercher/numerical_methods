#------------------------------------------
# Student: Andrew Kercher
# CSI 722: Project 1 for Fall 2011
# Last Updated: 12/15/2011
# Program: Makefile 
#------------------------------------------

#-----------------------
# VARIABLES: FORTRAN 90
#-----------------------

#---------------------------------
# External Modules for Fortran 90
#---------------------------------
# OBJ =  csi722_typeDef.o csi722_dataIO.o csi722_dataStruc.o csi722_sub.o csi722_project1.o 
MODF1 = csi722_typeDef.o
MODF2 = csi722_dataIO.o
MODF3 = csi722_dataStruc.o
MODF4 = csi722_elementOps.o
MODF5 = csi722_solveOps.o
SRCM1 = csi722_typeDef.f90
SRCM2 = csi722_dataIO.f90
SRCM3 = csi722_dataStruc.f90
SRCM4 = csi722_elementOps.f90
SRCM5 = csi722_solveOps.f90

#---------------------------------------------------
# Program and Output Files for csi721_project1.f90
#---------------------------------------------------
SOLF1 = data/pressure.zfem data/velocity.zfem
SOLF4 = data/euler_04.vtk
SOLF8 = data/euler_08.vtk
EXEF = csi722_project1.run
OBJF = csi722_project1.o 
SRCF = csi722_project1.f90

#--------------
# VARIABLES: C
#--------------
EXEC = csi722.run
OBJC = csi722_visual.o
SRCC = csi722_visual.c

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

# all: $(SOLF4)
# 	python remove_whitespace.py $?

04: #$(SOLF4)
	python remove_whitespace.py $(SOLF4)

08: #$(SOLF8)
	python remove_whitespace.py $(SOLF8)


#---------------------------------------------
# DEPENDENCIES: Compile needed *.f90 programs
#---------------------------------------------
$(MODF1): $(SRCM1)
	$(CMPF) $?

$(MODF2): $(SRCM2)
	$(CMPF) $?

$(MODF3): $(SRCM3)
	$(CMPF) $?

$(MODF4): $(SRCM4)
	$(CMPF) $?

$(MODF5): $(SRCM5)
	$(CMPF) $?

$(OBJF): $(SRCF) 
	$(CMPF) $?

$(EXEF): $(MODF1) $(MODF2) $(MODF3) $(MODF4) $(MODF5) $(OBJF)
	$(LNKF) -o $@ $(MODF1) $(MODF2) $(MODF3) $(MODF4) $(MODF5) $(OBJF) $(LIBF)

#--------------------------------------------------
# FIRST TARGET: Solve for scaler and vector fields
#--------------------------------------------------
# $(SOLF): $(EXEF)
# 	./$(EXEF)

#------------------------------------------------------------
# SECOND TARGET: Produce visualization using openGL and glut
#------------------------------------------------------------
$(EXEC): $(OBJC)
	$(LNKC) -o $@ $(OBJC) $(LIBC)

$(OBJC): $(SRCC)
	$(CMPC) $(SRCC)

#----------------------------------
# SECOND TARGET: Visualize Results
#----------------------------------
visual: $(EXEC) 
	./$? $(SOLF1)

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
