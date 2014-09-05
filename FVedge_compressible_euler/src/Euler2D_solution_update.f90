!!$ Program: Euler2D_solution_updates.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 02/09/2012
!!$ Description: Module includes subroutines to calculate 2D Euler fluxes.
!!$ Acknowledgment: The following algorithms and diagrams were taken from 
!!$                 or heavily based off of code originally written by 
!!$                 Dr. Katate Masatsuka (info[at]cfdbooks.com).


MODULE solution_update

private                        

public :: rk_update

CONTAINS

  !********************************************************************************
  !* -- Runge-Kutta explicit update --
  !*
  !* Input:           Crk = Runge-Kutta coefficient 
  !*                  dtg = global time step (not used if local time stepping)
  !*                   Cr = Courant number
  !*          nodes(:)res = residual at nodes
  !*
  !* Output:  ucons = updated conservative variables 
  !*          wprim = updated primitive variables 
  !*
  !********************************************************************************
  SUBROUTINE rk_update(Crk,dtg)

    USE constants
    USE mesh_data, only : nodes,ucons,wprim
    USE state_variables

    IMPLICIT NONE

    !Input variables:
    DOUBLE PRECISION,INTENT(IN) :: Crk,dtg

    !Local variables
    INTEGER :: i
    DOUBLE PRECISION,DIMENSION(nunks) :: du

    du(:) = 0.0d0
    do i = 1,nnode

       !calculate change in unknowns
       if(trim(step_type) == 'global')then
          !use global time-step
          du(:) = (Cr*Crk*dtg/nodes(i)%ml)*nodes(i)%res(:)
       elseif(trim(step_type) == 'local')then
          !local update
          du(:) = (Cr*Crk*nodes(i)%dt/nodes(i)%ml)*nodes(i)%res(:)
       endif

       !update solution
       ucons(i)%d  = ucons(i)%d + du(1) 
       ucons(i)%mx  = ucons(i)%mx + du(2) 
       ucons(i)%my  = ucons(i)%my + du(3) 
       ucons(i)%en  = ucons(i)%en + du(4) 

       !update primitive variables
       call cons2prim(i,i)         

    end do

  END SUBROUTINE rk_update

END MODULE solution_update
