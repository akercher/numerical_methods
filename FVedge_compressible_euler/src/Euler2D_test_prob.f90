!!$ Program: Euler2D_test_prob.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 02/09/2012
!!$ Description: Module includes initial solution of tests problems:
!!$              1. 2D shock diffraction
!!$              2. 2D naca0012 airfoil
!!$              2. 2D Kelvin Helmholtz instability

MODULE test_prob

private                        

public :: init_shock_diffraction
public :: init_airfoil
public :: init_kelvin_helmholtz
public :: init_shock_tube_2D

CONTAINS

!********************************************************************************
!********************************************************************************
!********************************************************************************
!* Initial state for 2D shock diffraction
!*
!********************************************************************************
  subroutine init_shock_diffraction

    USE constants
    USE mesh_data, only : wprim,ucons
    USE dataIO, only : printPrim,printCons
    USE state_variables, only : prim2cons

    IMPLICIT NONE

    INTEGER          :: i
    DOUBLE PRECISION :: M_shock, vx_shock, d0, vx0, vy0, pg0

    !initialize variables
    allocate(ucons(nnode),wprim(nnode))

    do i = 1,nnode

       ! Pre-shock state: uniform state; no disturbance has reahced yet.
       d0  = 1.0d0
       vx0 = 0.0d0
       vy0 = 0.0d0
       pg0 = 1.0d0/gamma

       ! Incoming shock speed

       M_shock = 5.090d0
       vx_shock = M_shock * sqrt(gamma*pg0/d0)

       ! Post-shock state: These values will be used in the inflow boundary condition.
       dinf = d0 * (gamma + 1.0d0)*M_shock**2/( (gamma - 1.0d0)*M_shock**2 + 2.0d0 )
       pginf =   pg0*(2.0d0*gamma*M_shock**2 - (gamma - 1.0d0) )/(gamma + 1.0d0)
       vxinf = (1.0d0 - d0/dinf)*vx_shock
       Minf  = vxinf / dsqrt(gamma*pginf/dinf)
       vyinf = 0.0d0
       
       ! Set the initial solution: set the pre-shock state inside the domain.
       wprim(i)%d  = d0
       wprim(i)%vx = vx0
       wprim(i)%vy = vy0
       wprim(i)%pg = pg0


       !calculate conservative variables
       call prim2cons(i,i)

    enddo
    
  end subroutine init_shock_diffraction
!--------------------------------------------------------------------------------

!********************************************************************************
!********************************************************************************
!********************************************************************************
!* Initial state for naca0012 airfoil
!*
!********************************************************************************

  subroutine init_airfoil

    USE constants
    USE mesh_data, only : wprim,ucons
    USE dataIO, only : printPrim,printCons
    USE state_variables, only : prim2cons

    IMPLICIT NONE

    INTEGER          :: i
!!$    DOUBLE PRECISION :: M_shock, vx_shock, d0, vx0, vy0, pg0

    !initialize variables
    allocate(ucons(nnode),wprim(nnode))

    !define far-field states
!!$    Minf = 0.40d0 !far-field Mach Number
    pginf = 1.0d0 
    csinf = 1.0d0 
    dinf = gamma*pginf/(csinf*csinf)
    vinf = Minf*csinf 
    vxinf = vinf*COS(alpha)
    vyinf = vinf*SIN(alpha)
    v2inf = vxinf*vxinf + vyinf*vyinf
    eninf = pginf/(gamma - 1.0d0) + dinf*v2inf/2.0d0    


    do i = 1,nnode
       
       !set all values to far-field values
       wprim(i)%d  = dinf
       wprim(i)%vx = vxinf
       wprim(i)%vy = vyinf
       wprim(i)%pg = pginf

       !calculate conservative variables
       call prim2cons(i,i)

    enddo
    
  end subroutine init_airfoil

!--------------------------------------------------------------------------------

!********************************************************************************
!********************************************************************************
!********************************************************************************
!* Initial conditions for Kelvin-Helmholtz instability
!*
!********************************************************************************

  subroutine init_kelvin_helmholtz

    USE constants
    USE mesh_data, only : nodes,wprim,ucons
    USE dataIO, only : printPrim,printCons
    USE state_variables, only : prim2cons

    IMPLICIT NONE

    INTEGER          :: i
    DOUBLE PRECISION :: M1,M2   !mach numbers
    DOUBLE PRECISION :: d1,d2   !densities
    DOUBLE PRECISION :: v0,vx1,vx2,vy   !velocitys
    DOUBLE PRECISION :: pg                !pressure
    DOUBLE PRECISION :: cs1,cs2           !speed of sound
    DOUBLE PRECISION :: ap                !preturbation amplitude
    DOUBLE PRECISION :: y1,y2             !interface locations
    DOUBLE PRECISION :: xi                !x component of node i
!!$    DOUBLE PRECISION :: nx
!!$    DOUBLE PRECISION,allocatable(:) :: rannum

    !set state variables and perturbation amplitude
    d1 = 1.0d0
    d2 = 2.0d0
    v0 = 0.50d0
    pg = 2.50d0
    cs1 = dsqrt(gamma*pg/d1)
    cs2 = dsqrt(gamma*pg/d2)
    M1 = v0/cs1
    M2 = v0/cs2

    ap = 0.010d0

    !initialize variables
    allocate(ucons(nnode),wprim(nnode))

    !set interface locations
    y1 = Ly/4.0d0
    y2 = 3.0d0*Ly/4.0d0

    do i = 1,nnode
       xi = nodes(i)%x
       vy = -ap*DSIN(2.0d0*pi*xi/Lx)
       if(nodes(i)%y < y1 .or. nodes(i)%y > y2)then
          !set to state 1
          wprim(i)%d  = d1
          wprim(i)%vx = 2.0*(v0 + vy) 
          wprim(i)%vy = vy
          wprim(i)%pg = pg
       else
          !set to state 2
          wprim(i)%d  = d2
          wprim(i)%vx = -2.0*(v0 + vy) 
          wprim(i)%vy = vy
          wprim(i)%pg = pg          
       endif

       !calculate conservative variables
       call prim2cons(i,i)

    enddo
    
  end subroutine init_kelvin_helmholtz
 !--------------------------------------------------------------------------------

!********************************************************************************
!********************************************************************************
!********************************************************************************
!* Initial conditions for 2D shock tube
!*
!********************************************************************************

  subroutine init_shock_tube_2D

    USE constants
    USE mesh_data, only : nodes,wprim,ucons
    USE dataIO, only : printPrim,printCons
    USE state_variables, only : prim2cons

    IMPLICIT NONE

    INTEGER          :: i
    DOUBLE PRECISION :: Ml,Mr   !mach numbers
    DOUBLE PRECISION :: dl,dr   !densities
    DOUBLE PRECISION :: vnl,vnr !velocitys normal to shock
    DOUBLE PRECISION :: vxl,vyl !left-velocitys
    DOUBLE PRECISION :: vxr,vyr !right-velocitys
    DOUBLE PRECISION :: pgl,pgr !pressures
    DOUBLE PRECISION :: csl,csr !speed of sound
    DOUBLE PRECISION :: xs      !interface location
    DOUBLE PRECISION :: xi      !x component of node i
!!$    DOUBLE PRECISION :: nx
!!$    DOUBLE PRECISION,allocatable(:) :: rannum

    !set state variables
    dl  = 1.0d0
    vnl = 0.750d0
    vxl = vnl*cos(alpha)
    vyl = vnl*sin(alpha)
    pgl = 1.0d0
    csl = sqrt(gamma*pgl/dl)

    dr  = 0.1250d0
    vnr = 0.0
    vxr = vnr*cos(alpha)
    vyr = vnr*sin(alpha)
    pgr = 0.10d0
    csr = sqrt(gamma*pgr/dr)

    Ml = vnl/csl
    Mr = vnr/csr

    !initialize variables
    allocate(ucons(nnode),wprim(nnode))

    !set interface locations
    xi = Lx/2.0d0

    do i = 1,nnode

       if(nodes(i)%x .le. xi )then
          !set to state 1
          wprim(i)%d  = dl
          wprim(i)%vx = vxl
          wprim(i)%vy = vyl
          wprim(i)%pg = pgl
       else
          !set to state 2
          wprim(i)%d  = dr
          wprim(i)%vx = vxr
          wprim(i)%vy = vyr
          wprim(i)%pg = pgr          
       endif

       !calculate conservative variables
       call prim2cons(i,i)

!!$       call printPrim(i,i)

    enddo
    
  end subroutine init_shock_tube_2D
 !--------------------------------------------------------------------------------


END MODULE test_prob
