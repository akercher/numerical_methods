#------------------------------------------
# Author:       Andrew Kercher
# Program:      Makefile
# Last Updated: 02/09/2012 
# Description:  Makefile for 2D Euler Eqs.
#------------------------------------------

#-----------------------
# VARIABLES: FORTRAN 90
#-----------------------

#-------------------------------------
# Source Files for Fortran 90 Modules
#-------------------------------------
SRCS = src/Euler2D_defs.f90              \
       src/Euler2D_dataIO.f90            \
       src/Euler2D_initCalc.f90          \
       src/Euler2D_initMesh.f90          \
       src/Euler2D_timeCLF.f90           \
       src/Euler2D_state_variables.f90   \
       src/Euler2D_fluxes.f90            \
       src/Euler2D_BCs.o                 \
       src/Euler2D_solution_update.f90   \
       src/Euler2D_solver.f90            \
       src/FV_mesh_check.f90             \
       src/Euler2D_test_prob.f90         \

#------------------------
# Objests for Fortran 90
#------------------------
OBJS = Euler2D_defs.o              \
       Euler2D_dataIO.o            \
       Euler2D_initCalc.o          \
       Euler2D_initMesh.o          \
       Euler2D_timeCLF.o           \
       Euler2D_state_variables.o   \
       Euler2D_fluxes.o            \
       Euler2D_BCs.o               \
       Euler2D_solution_update.o   \
       Euler2D_solver.o            \
       FV_mesh_check.o             \
       Euler2D_test_prob.o         \

#------------------------------------------------
# Program and Output Files for fe_problem.f90
#------------------------------------------------
SOLF1 = data/pressure.zfem data/velocity.zfem
SOLF4 = data/euler_04.vtk
SOLF8 = data/euler_08.vtk

#------------------------
# Problem Specific Files
#------------------------
PROB = Euler2D_main
PROBEXE = $(PROB).run
PROBOBJ = $(PROB).o
PROBSRC = src/$(PROB).f90

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
example: $(PROBEXE)
	./$?


#---------------------------------------------
# DEPENDENCIES: Compile needed *.f90 programs
#---------------------------------------------
%.mod: %.o
%.o: ./src/%.f90
	$(CMPF) $<

$(PROBOBJ): $(OBJS)
Euler2D_dataIO.o          : Euler2D_defs.o
Euler2D_timeCLF.o         : Euler2D_defs.o
Euler2D_state_variables.o : Euler2D_defs.o
Euler2D_fluxes.o          : Euler2D_defs.o Euler2D_state_variables.o
Euler2D_BCs.o             : Euler2D_defs.o Euler2D_state_variables.o Euler2D_fluxes.o 
Euler2D_solution_update.o : Euler2D_defs.o Euler2D_state_variables.o
Euler2D_solver.o          : Euler2D_defs.o Euler2D_dataIO.o Euler2D_timeCLF.o Euler2D_state_variables.o Euler2D_fluxes.o Euler2D_BCs.o Euler2D_solution_update.o
Euler2D_initCalc.o        : Euler2D_defs.o Euler2D_dataIO.o 
Euler2D_initMesh.o        : Euler2D_defs.o Euler2D_solver.o 
FV_mesh_check.o           : Euler2D_defs.o Euler2D_solver.o 
Euler2D_test_prob.o       : Euler2D_defs.o Euler2D_dataIO.o Euler2D_state_variables.o 

$(PROBEXE): $(OBJS) $(PROBOBJ)
	$(LNKF) -o $@ $? $(LIBF) 

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
	rm -f src/*.o
	rm -f src/*.mod
	rm -f src/*~
	rm -f *.o
	rm -f *.mod
	rm -f *.run
	rm -f *~