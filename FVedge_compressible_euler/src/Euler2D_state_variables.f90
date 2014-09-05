!!$ Program: Euler2D_state_variables.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 02/09/2012
!!$ Description: Module includes subroutines to convert conservative variables to
!!$              primative variables and vise-versa

MODULE state_variables

private                        

public :: cons2prim
public :: prim2cons


CONTAINS


!********************************************************************************
!* Compute primative variables from conservative variables
!*
!* ------------------------------------------------------------------------------
!*  Input:  ucons = conservative variables (d, d*vx,  d*vy, d*E)
!* Output:  wprim =    primitive variables (d,    vx,   vy,   p)
!* ------------------------------------------------------------------------------
!* 
!********************************************************************************
  subroutine cons2prim(istart,iend)

    USE constants, only : half,nunks,gamma,pgmin
    USE mesh_data, only : ucons,wprim

    IMPLICIT NONE
    
    INTEGER,INTENT(IN)                :: istart,iend !starting and ending node
    INTEGER                           :: i,j,k
    DOUBLE PRECISION                  :: KE !kinetic energy

    do i=istart,iend
       KE          = half*(ucons(i)%mx*ucons(i)%mx + ucons(i)%my*ucons(i)%my)/ucons(i)%d
       wprim(i)%d  = ucons(i)%d
       wprim(i)%vx = ucons(i)%mx/ucons(i)%d
       wprim(i)%vy = ucons(i)%my/ucons(i)%d
       wprim(i)%pg = max(pgmin,(gamma-1.0d0)*(ucons(i)%en - KE))
    enddo

  end subroutine cons2prim

!********************************************************************************
!* Compute conservative variables from primative variables
!*
!* ------------------------------------------------------------------------------
!*  Input:  wprim =    primitive variables (d,    vx,   vy,   p)
!* Output:  ucons = conservative variables (d, d*vx,  d*vy, d*E)
!* ------------------------------------------------------------------------------
!* 
!********************************************************************************
  subroutine prim2cons(istart,iend)

    USE constants, only : half,nunks,gamma,pgmin
    USE mesh_data, only : ucons,wprim

    IMPLICIT NONE
    
    INTEGER,INTENT(IN)                :: istart,iend !starting and ending node
    INTEGER                           :: i,j,k
    DOUBLE PRECISION                  :: KE !kinetic energy

    do i=istart,iend
       KE          = half*(wprim(i)%vx*wprim(i)%vx + wprim(i)%vy*wprim(i)%vy)
       ucons(i)%d  = wprim(i)%d
       ucons(i)%mx = wprim(i)%d*wprim(i)%vx
       ucons(i)%my = wprim(i)%d*wprim(i)%vy
       ucons(i)%en = wprim(i)%pg/(gamma - 1.0d0) + wprim(i)%d*KE
    enddo

  end subroutine prim2cons
  !--------------------------------------------------------------------------------

END MODULE state_variables
