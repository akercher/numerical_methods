!!$ Program: Euler2D_solver.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 02/09/2012
!!$ Description: Module includes subroutines for 2D Euler solvers.
!!$ Acknowledgment: The following algorithms and diagrams were taken from 
!!$                 or heavily based off of code originally written by 
!!$                 Dr. Katate Masatsuka (info[at]cfdbooks.com).


MODULE solver

private                        

public :: FVedge_solver
public :: lsq_inv_2x2
public :: lsq_grads_2x2

CONTAINS

!********************************************************************************
!********************************************************************************
!********************************************************************************
!* Euler solver: Node-Centered Finite-Volume Method (Edge-Based)
!*
!* - Node-centered finite-volume method for unstructured grids(quad/tri/mixed)
!* - Roe flux with an entropy fix and Rotated-RHLL flux
!* - Reconstruction by unweighted least-squares method (2x2 system for gradients)
!* - Van Albada slope limiter to the primitive variable gradients
!* - 2-Stage Runge-Kutta time-stepping
!*
!********************************************************************************
  subroutine FVedge_solver

    use constants
    use data_types, only : cons
    use dataIO, only : printCons,printPrim
    use mesh_data, only : nodes,bouns,ucons,wprim,gradw
    use timeCLF, only : local_time_step
    use BCs, only : wall_normal_mass_flux
    use solution_update

    IMPLICIT NONE

    INTEGER                                     :: i,j,k
    INTEGER                                     :: kstep    !number of steps
    INTEGER                                     :: nstep    !number of iterations
    DOUBLE PRECISION                            :: dtg      !global minimum time-step
    DOUBLE PRECISION                            :: time     !current time
    DOUBLE PRECISION,DIMENSION(nunks,3)         :: res_norm !Residual norms(L1,L2,Linf)
    TYPE(CONS),DIMENSION(:),POINTER             :: u0       !solution stored for rk update
    DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: du       !change in cons. variables



    !allocate the temporary solution array needed for the Runge-Kutta method
    !and change in conservative variables
    allocate(u0(nnode),du(nunks))

    ! These parameters are set in main.f90. Here just print them on display.
    write(*,*)
    write(*,*) "Calling the Euler solver..."
    write(*,*)
    write(*,*) "                     Minf = ", Minf
    write(*,*) "           Courant number = ", Cr
    write(*,*) "               final time = ", tf
    write(*,*) "     max. number of steps = ", max_steps
    write(*,*) "                flux_type = ", trim(flux_type)
    write(*,*) "             limiter_type = ", trim(limiter_type)
    write(*,*)

    !--------------------------------------------------------------------------------
    ! First, make sure that normal mass flux is zero at all solid boundary nodes.
    ! NOTE: Necessary because initial solution may generate the normal component.
    !--------------------------------------------------------------------------------

    call wall_normal_mass_flux

    !--------------------------------------------------------------------------------
    ! Time-stepping toward the final time
    !--------------------------------------------------------------------------------
    time = 0.0d0

    !initialize gradients for reconstruction
    allocate(gradw(nnode)) 
    steps : do kstep = 1,max_steps

       !------------------------------------------------------
       ! Two-stage Runge-Kutta scheme: u^n is saved as u0(:,:)
       !  1. u^*     = u^n - (dt/vol)*res(u^n)
       !  2. u^{n+1} = 1/2*u^n + 1/2*[u^* - (dt/vol)*res(u^*)]
       !------------------------------------------------------

       !-----------------------------
       !- 1st Stage of Runge-Kutta:
       !-----------------------------

       !compute res(u^n)
       call fvedge_residual

       !compute residual norms (divided residual)
       call residual_norm(res_norm)

       !print the residual norm.
       if (kstep==1) write(*,'(a88)') "Density    X-momentum  Y-momentum   Energy"
       if (mod(kstep,10)==0) write(*,'(a2,es9.3,a9,i10,a12,4es12.2)') &
            "t=", time, "steps=", kstep, " L1(res)=",res_norm(:,1)

       if(trim(step_type) == 'global')then
          !   Stop if the final time is reached.
          if (time == tf) exit steps
       end if

       !save current solution for later stage rk step
       do i = 1, nnode
          u0(i)%d  = ucons(i)%d
          u0(i)%mx = ucons(i)%mx
          u0(i)%my = ucons(i)%my
          u0(i)%en = ucons(i)%en
       end do
       
       !Compute the time step (local and global)
       call local_time_step(dtg)
       
       !Adjust dt so as to finish exactly at the final time
       if (time + dtg > tf) dtg = (tf - time)

       !Update the solution
       !1st Stage => u^* = u^n - dt/dx*Res(u^n)
       call rk_update(1.0d0,dtg)

       !-----------------------------
       !- 2nd Stage of Runge-Kutta:
       !-----------------------------

       !compute Res(u^*)
       call fvedge_residual
       
       !compute 1/2*(u^n + u^*)
       do i = 1,nnode
          ucons(i)%d  = half*(ucons(i)%d + u0(i)%d)
          ucons(i)%mx = half*(ucons(i)%mx + u0(i)%mx)
          ucons(i)%my = half*(ucons(i)%my + u0(i)%my)
          ucons(i)%en = half*(ucons(i)%en + u0(i)%en)
       end do

       !2nd Stage => u^{n+1} = 1/2*(u^n + u^*) - 1/2*dt/dx*Res(u^*)
       call rk_update(half,dtg)

       time = time + dtg

    end do steps

    !--------------------------------------------------------------------------------
    ! End of Time-stepping to the final time.
    !--------------------------------------------------------------------------------
    
    if (kstep == max_steps .and. (time /= tf)) then
       write(*,*) " Final time not reached... Sorry..."
       write(*,*) "   Max time step reached... time_step_max=", max_steps
       write(*,*) "   Increase time_step_max, and try again."
       write(*,*) " ... Stop"
       stop
    endif

    write(*,*) " Congratulations. Final time is reached."
    write(*,*) " Final step number     =", kstep
    write(*,*) " Final time           = ", time
    write(*,*) " Final time requested = ", tf
    write(*,*)
    write(*,*) "Finished the Euler solver... Bye!"
    
  end subroutine FVedge_solver

  !********************************************************************************
  !* This subroutine computes the residual for a node-centered finite-volume method
  !*
  !* ------------------------------------------------------------------------------
  !*  Input: the current solution
  !*
  !* Output: node(:)%res = the residual computed by the current solution.
  !* ------------------------------------------------------------------------------
  !*
  !* Note: dU/dt + dF/dx + dG/dy = 0. Residuals are first computed as
  !*       the integral of (dF/dx + dG/dy), and at the end negative sign is added
  !*       so that we have dU/dt = Res at every node.
  !* 
  !********************************************************************************
  subroutine fvedge_residual
    
    use constants
    use mesh_data, only : nodes,edges,bouns,elems,ucons,wprim,gradw
    use fluxes,    only : fluxRoe2D,fluxRHLL2D,fluxHLLC2D
    use state_variables
    use BCs

    IMPLICIT NONE
    
    INTEGER                            :: i,j,k
    INTEGER                            :: nel, ner   !left/right nodes of each edge
    INTEGER                            :: ni,nj   !left/right nodes of each edge
    INTEGER                            :: belem   !element adjacent to boundary face
    INTEGER                            :: nfl,nfr !left/right nodes of boundary face
    INTEGER                            :: ix=1,iy=2
    DOUBLE PRECISION                   :: se        !magnitude of the edge vector
    DOUBLE PRECISION                   :: sa        !magnitude of the area-normal vector
    DOUBLE PRECISION                   :: bfn_mag        !magnitude of boundary face normal
    DOUBLE PRECISION                   :: nwsd     !normal maximum wave speed
    DOUBLE PRECISION                   :: mnorm    !normal component of the momentum
    DOUBLE PRECISION, DIMENSION(nunks) :: Fedge    !flux along the edge
    DOUBLE PRECISION, DIMENSION(nunks) :: wl,wr      !left/right primative states
    DOUBLE PRECISION, DIMENSION(nunks) :: wva
    DOUBLE PRECISION, DIMENSION(nunks) :: dwl,dwr       !slope at left/right nodes
    DOUBLE PRECISION, DIMENSION(nunks) :: dwm, dwp, dwe !re-defined slopes to be limited
    DOUBLE PRECISION, DIMENSION(2)     :: sev            !unit edge vector
    DOUBLE PRECISION, DIMENSION(2)     :: sav            !unit directed area vector
    DOUBLE PRECISION, DIMENSION(2)     :: bfn            !unit boundary face normal
    DOUBLE PRECISION, DIMENSION(2)     :: nbd            !unit boundary node normal
    DOUBLE PRECISION, DIMENSION(nunks) :: bfl, bfr !boundary flux at left/right nodes

    !-------------------------------------------------------------------------
    ! Gradient Reconstruction for second-order accuracy
    
    !  Initialization
    
    nodes0 : do i = 1,nnode
       gradw(i)%dx(:) = 0.0d0
       gradw(i)%dy(:) = 0.0d0
    end do nodes0
    
    !  Perform gradient reconstruction
    nodes1 : do i = 1,nnode
       call lsq_grads_2x2(i)
    end do nodes1
    
    
    !-------------------------------------------------------------------------
    ! Residual computation: interior fluxes
    
    !initialization
    do i = 1,nnode
       nodes(i)%res(:) = 0.0d0
       nodes(i)%nwsd   = 0.0d0
    enddo
    
    ! Flux computation across internal edges (to be accumulated in res(:))
    !
    !      nj              1. Extrapolate the solutions to the edge-midpoint
    !       o                 from the nodes, ni and nj.
    !        \   face      2. Compute the numerical flux
    !         \ -------c2  3. Add it to the residual for ni, and subtract it from
    !        / \              the residual for nj.
    !   face/   \ edge
    !      /     o         Directed area is the sum of the left and the right faces.
    !    c1       ni       Left/right face is defined by the edge-midpoint and
    !                      the centroid of the left/right element.
    !                      Directed area is positive in ni -> nj
    !
    ! (c1, c2: element centroids)
    !
    !--------------------------------------------------------------------------------
    edge_loop : do i = 1,nedge
       !--------------------------------------------------------------------------------
       
       ! Left and right nodes of the i-th edge
       
       nel = edges(i)%ni  ! left node of the edge
       ner = edges(i)%nj  ! right node of the edge
       sav = edges(i)%sav ! This is the directed face vector (unit vector)
       sa  = edges(i)%sa  ! Magnitude of the directed area vector
       sev = edges(i)%sev  ! This is the vector along the edge (unit vector)
       se  = edges(i)%se   ! Magnitude of the edge vector (Length of the edge)

       ! Solution gradient projected along the edge
       !
       !  NOTE: The gradient is multiplied by the distance.
       !        So, it is equivalent to the solution difference.
       
       dwl(:) = half*se*(gradw(nel)%dx(:)*sev(1) + gradw(nel)%dy(:)*sev(2))
       dwr(:) = half*se*(gradw(ner)%dx(:)*sev(1) + gradw(ner)%dy(:)*sev(2))
       
       !(1) No limiter (good for smooth solutions)
       limiter : if (trim(limiter_type) == "none") then
          
          !      Simple linear extrapolation
          wl(1) = wprim(nel)%d + dwl(1) 
          wl(2) = wprim(nel)%vx + dwl(2) 
          wl(3) = wprim(nel)%vy + dwl(3) 
          wl(4) = wprim(nel)%pg + dwl(4)
          
          wr(1) = wprim(ner)%d - dwr(1) 
          wr(2) = wprim(ner)%vx - dwr(2) 
          wr(3) = wprim(ner)%vy - dwr(3) 
          wr(4) = wprim(ner)%pg - dwr(4) 
          
          
       !(2) UMUSCL-type limiters: simple 1D limiting.   
       elseif (trim(limiter_type) == "vanalbada") then
          
          !In 1D: dwp = w_{j+1}-w_j, dwm = w_j-w_{j-1} => limited_slope = limiter(dwm,dwp)
          !       We can do the same in 2D as follows.
          !In 2D:    dwp = w_{neighbor}-w_j, dwm = 2*(grad_w_j*edge)-dwp
          !        => limited_slope = limiter(dwm,dwp)
          !NOTE: On a regular grid, grad_w_j*edge will be the central-difference,
          !      so that the average (dwm+dwp)/2 will be the central-difference just like in 1D.
          
          !edge derivative
          dwe(1) = half*(wprim(ner)%d - wprim(nel)%d)
          dwe(2) = half*(wprim(ner)%vx - wprim(nel)%vx)
          dwe(3) = half*(wprim(ner)%vy - wprim(nel)%vy)
          dwe(4) = half*(wprim(ner)%pg - wprim(nel)%pg)
          
          !left face value (wl) with the Van Albada limiter
          dwm(:)  = 2.0d0*dwl(:) - dwe(:)
          dwp(:)  = dwe(:)
          wva(:)  = va_limiter(dwm,dwp,se)
          
          wl(1) = wprim(nel)%d + wva(1)
          wl(2) = wprim(nel)%vx + wva(2)
          wl(3) = wprim(nel)%vy + wva(3)
          wl(4) = wprim(nel)%pg + wva(4)
          
          !right face value (wr) with the Van Albada limiter
          dwm(:)  = -(2.0d0*dwr(:) - dwe(:))
          dwp(:)  = -dwe(:)
          wva(:)  = va_limiter(dwm,dwp,se)
          
          wr(1) = wprim(ner)%d + wva(1)
          wr(2) = wprim(ner)%vx + wva(2)
          wr(3) = wprim(ner)%vy + wva(3)
          wr(4) = wprim(ner)%pg + wva(4)
          
       !(3) No other limiters available.
          
       else
         
          write(*,*) " Invalid input for limiter_type = ", trim(limiter_type)
          write(*,*) " Choose none or vanalbada, and try again."
          write(*,*) " ... Stop."
          stop
          
       endif limiter

       !initialize flux
       Fedge(:) = 0.0d0
       !calculate flux for given left/right states.      
       !(1) Roe flux (carbuncle is expected for strong shocks)
       if     (trim(flux_type)=="roe") then
          
          call fluxRoe2D(sav,wl,wr,Fedge,nwsd)
          
       !(2) Rotated-RHLL flux (no carbuncle is expected for strong shocks)
       elseif (trim(flux_type)=="rhll") then

          call fluxRHLL2D(sav,wl,wr,Fedge,nwsd)

       !(3) HLLC flux (resolves contact wave)
       elseif (trim(flux_type)=="hllc") then

          call fluxHLLC2D(sav,wl,wr,Fedge,nwsd)
          
       else
          
          write(*,*) " Invalid input for flux_type = ", trim(flux_type)
          write(*,*) " Choose roe, rhll or hllc, and try again."
          write(*,*) " ... Stop."
          stop
          
       endif
       
       !Add the flux multiplied by the magnitude of the directed area vector to left node,
       !and accumulate the max wave speed quantity for use in the time step calculation.

       nodes(nel)%res(:) = nodes(nel)%res(:)  +  Fedge(:)*sa
       nodes(nel)%nwsd   = nodes(nel)%nwsd  + nwsd*sa

       !Subtract the flux multiplied by the magnitude of the directed area vector from right node,
       !and accumulate the max wave speed quantity for use in the time step calculation.
       !
       !NOTE: Subtract because the outward face normal is -sa for the right node.
       
       nodes(ner)%res(:) = nodes(ner)%res(:)  -  Fedge(:)*sa
       nodes(ner)%nwsd   = nodes(ner)%nwsd  + nwsd*sa

       
       !--------------------------------------------------------------------------------
       
    end do edge_loop
    !--------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------
    !apply boundary conditions
    call FVedge_BCs
    !--------------------------------------------------------------------------------    
    
    !Switch the residual sign.   
    do i = 1,nnode   
       nodes(i)%res(:) = - nodes(i)%res(:)
    end do
    
  end subroutine fvedge_residual
  !--------------------------------------------------------------------------------
  
  !********************************************************************************
  !* -- vanAlbada Slope Limiter Function--
  !*
  !* 'A comparative study of computational methods in cosmic gas dynamics', 
  !* Van Albada, G D, B. Van Leer and W. W. Roberts, Astronomy and Astrophysics,
  !* 108, p76, 1982
  !*
  !* ------------------------------------------------------------------------------
  !*  Input:   da, db     : two differences
  !*
  !* Output:   va_limiter : limited difference
  !* ------------------------------------------------------------------------------
  !*
  !********************************************************************************
  pure function va_limiter(da,db,h)
    
    use constants   , only : half
    
    IMPLICIT NONE 
    
    DOUBLE PRECISION,               INTENT(in) :: h
    DOUBLE PRECISION, DIMENSION(4), INTENT(in) :: da, db
    DOUBLE PRECISION, DIMENSION(4)             :: va_limiter
    DOUBLE PRECISION :: eps2
    
    continue
    
    eps2 = (0.30d0*h)**3
    
    va_limiter = half*( sign(1.0d0,da*db) + 1.0d0 ) * &
         ( (db**2 + eps2)*da + (da**2 + eps2)*db )/(da**2 + db**2 + 2.0d0*eps2)
    
    
  end function va_limiter
  !--------------------------------------------------------------------------------

  
  !********************************************************************************
  !********************************************************************************
  !********************************************************************************
  !* Prepararion for Tangency condition (slip wall):
  !*
  !* Eliminate normal mass flux component at all solid-boundary nodes at the
  !* beginning. The normal component will never be changed in the solver: the
  !* residuals will be constrained to have zero normal component.
  !*
  !********************************************************************************
!!$  subroutine wall_normal_mass_flux
!!$    
!!$    use constants, only : nunks,nbcs,prob_type
!!$    use mesh_data, only : nodes,bouns,ucons,wprim
!!$    use state_variables
!!$
!!$    IMPLICIT NONE
!!$
!!$    INTEGER                       :: i,j,inode
!!$    DOUBLE PRECISION              :: nmf     !normal mass flux
!!$    DOUBLE PRECISION,DIMENSION(2) :: nbd     !boundary normal
!!$
!!$    bc_loop : do i = 1, nbcs
!!$
!!$       only_slip_wall : if(trim(bouns(i)%bname) == "slip_wall") then
!!$
!!$          write(*,*) " Eliminating the normal momentum on slip wall boundary ", i
!!$
!!$          bnode_slip_wall : do j = 1,bouns(i)%nbnode
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$! THIS IS A SPECIAL TREATMENT FOR SHOCK DIFFRACTION PROBLEM.
!!$!
!!$! NOTE: This is a corner point between the inflow boundary and
!!$!       the lower-left wall. Enforce zero y-momentum, which is
!!$!       not ensured by the standard BCs.
!!$!       This special treatment is necessary because the domain
!!$!       is rectangular (the left boundary is a straight ine) and
!!$!       the midpoint node on the left boundary is actually a corner.
!!$!
!!$!       Our computational domain:
!!$!
!!$!                 ---------------
!!$!          Inflow |             |
!!$!                 |             |  o: Corner node
!!$!          .......o             |
!!$!            Wall |             |  This node is a corner.
!!$!                 |             |
!!$!                 ---------------
!!$!
!!$!       This is to simulate the actual domain shown below:
!!$!      
!!$!         -----------------------
!!$! Inflow  |                     |
!!$!         |                     |  o: Corner node
!!$!         --------o             |
!!$!            Wall |             |
!!$!                 |             |
!!$!                 ---------------
!!$!      In effect, we're simulating this flow by a simplified
!!$!      rectangular domain (easier to generate the grid).
!!$!      So, an appropriate slip BC at the corner node needs to be applied,
!!$!      which is "zero y-momentum", and that's all.
!!$!
!!$             if(trim(prob_type) == 'shock_diff')then
!!$                if (i==2 .and. j==1) then
!!$                   inode = bouns(i)%bnode(j)
!!$                   ucons(inode)%my = 0.0d0     !set my to zero
!!$                   call cons2prim(inode,inode) !update primitive variables
!!$                   cycle bnode_slip_wall       !go to the next.
!!$                endif
!!$             endif
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$             inode = bouns(i)%bnode(j)
!!$             nbd(1) = bouns(i)%bnx(j)
!!$             nbd(2) = bouns(i)%bny(j)
!!$             
!!$             !inner product of momentum and normal
!!$             nmf = ucons(inode)%mx*nbd(1) + ucons(inode)%my*nbd(2)
!!$
!!$             !eliminate normal mass flux from boundary
!!$             ucons(inode)%mx = ucons(inode)%mx - nmf*nbd(1)
!!$             ucons(inode)%my = ucons(inode)%my - nmf*nbd(2)
!!$             
!!$             call cons2prim(inode,inode) !Update primitive variables
!!$
!!$          end do bnode_slip_wall
!!$
!!$          write(*,*) " Finished eliminating the normal momentum on slip wall boundary ", i
!!$          write(*,*)
!!$          
!!$       endif only_slip_wall
!!$       
!!$    end do bc_loop
!!$
!!$  end subroutine wall_normal_mass_flux
!--------------------------------------------------------------------------------


!********************************************************************************
!* Given the primitive variables at nodes {k} around node i, this computes the 
!* gradients (wx,wy) at node j by the unweighted least-squares method.
!*
!* ------------------------------------------------------------------------------
!*  Input:  wprim(inode) = current primitive variables at node 'inode'
!*
!* Output:  gradw(inode) = gradients of the primitive variables
!* ------------------------------------------------------------------------------
!*  
!*  1. At node i, compute a vector b = \sum_k [ (xk-xi)*(uk-ui), (yk-yi)*(uk-ui) ]
!*  2. Compute the gradient by multiplying b by the inverse LSQ matrix that has
!*     been pre-computed by the subroutine lsq01_matrix() in the main:
!*             ux = inverse(1,1)*b(1) + inverse(1,2)*b(2)
!*             uy = inverse(2,1)*b(1) + inverse(2,2)*b(2)
!*
!********************************************************************************
  subroutine lsq_grads_2x2(inode)

    USE constants
    USE mesh_data, only : nodes,wprim,gradw

    IMPLICIT NONE
    
    INTEGER,INTENT(IN)                :: inode
    INTEGER                           :: k,in,isuno
    INTEGER                           :: ix,iy
    DOUBLE PRECISION                  :: dvar,dx,dy
    DOUBLE PRECISION,DIMENSION(2)     :: B
    DOUBLE PRECISION,DIMENSION(nunks) :: wi,wsuno
    
    ix=1 
    iy=2

!!$print*,'inode =',inode
!!$print*,nodes(inode)%lsq_inv(ix,1),nodes(inode)%lsq_inv(ix,2)
!!$print*,nodes(inode)%lsq_inv(iy,1),nodes(inode)%lsq_inv(iy,2)    

    !primative variables at current node
    wi(:) = 0.0d0
    wi(1) = wprim(inode)%d
    wi(2) = wprim(inode)%vx
    wi(3) = wprim(inode)%vy
    wi(4) = wprim(inode)%pg

    !primative variables at surrounding node
    wsuno(:) = 0.0d0

    ! Loop over variables
    do k = 1, 4

       B(:) = 0.0d0

       !Loop over surrounding elements
       do in = 1, nodes(inode)%nsuno
          isuno = nodes(inode)%suno(in)
       
          wsuno(1) = wprim(isuno)%d
          wsuno(2) = wprim(isuno)%vx
          wsuno(3) = wprim(isuno)%vy
          wsuno(4) = wprim(isuno)%pg

          dvar = wsuno(k) - wi(k) !primitive variables

          dx = nodes(isuno)%x    -  nodes(inode)%x
          dy = nodes(isuno)%y    -  nodes(inode)%y

          B(1) = B(1) + dvar*dx
          B(2) = B(2) + dvar*dy

!!$print*,k,isuno,inode
!!$print*,k,wsuno(k),wi(k)
          
       enddo

!!$print*,k,B(1),B(2)

       !  Multiply the inverse LSQ matrix to get the gradients
       gradw(inode)%dx(k) = nodes(inode)%lsq_inv(ix,1)*B(1) &
            + nodes(inode)%lsq_inv(ix,2)*B(2)
       
       gradw(inode)%dy(k) = nodes(inode)%lsq_inv(iy,1)*B(1) &
            + nodes(inode)%lsq_inv(iy,2)*B(2)

    enddo

  END SUBROUTINE lsq_grads_2x2

!--------------------------------------------------------------------------------
!********************************************************************************
!* --- Inverse Matrix for 2x2 Least-Squares Gradient Reconstruction ---
!*
!* Construct a matrix for the linear least-squares(LSQ) gradient reconstruction.
!* (unweighted LSQ; more accurate than weighted ones to my knowledge.)
!*
!* Note: it requires at least 2 non-colinear neighbors.
!*
!* Example: Consider constructing (ux,uy) at i with the following stencil.
!*
!*      3 o     o 2
!*         \   / 
!*          \ /
!*         i *-----o 1
!*          /|
!*         / |
!*        /  o 5      *: node in interest (i)
!*       o 4          o: neighbors (k = 1,2,3,4,5)
!*
!*  5 equations:
!*    (x1-xi)*ux + (y1-yi)*uy = (u1-ui)
!*    (x2-xi)*ux + (y2-yi)*uy = (u2-ui)
!*    (x3-xi)*ux + (y3-yi)*uy = (u3-ui)
!*    (x4-xi)*ux + (y4-yi)*uy = (u4-ui)
!*    (x5-xi)*ux + (y5-yi)*uy = (u5-ui)
!*
!*  This system is written in the matrix form:
!*
!*        A*x = b,  x=(ux,uy), A=5x2 matrix, b=5x1 matrix
!*
!*  The least-squares problem is
!*
!*      A^T*A*x = A^T*b, (T denotes the transpose: A^T=2x5 matrix)
!*  
!*  which is
!*
!*  [sum_k (xk-xi)^2]*ux       + [sum_k (xk-xi)*(yk-yi)]*uy = [sum_k (uk-ui)*(xk-xi)]
!*  [sum_k (xk-xi)*(yk-yi)]*ux + [sum_k (yk-yi)]*uy         = [sum_k (uk-ui)*(yk-yi)]
!*
!* This subroutine computes the inverse of (A^T*A) at every node (which depends
!* only on the grid geometry), so that the gradient at a node can be computed
!* by a matrix-vector multiplication, i.e., (A^T*A)^{-1}*(A^T*b), 
!* (only A^T*b needs to be re-computed).
!*
!* ------------------------------------------------------------------------------
!*  Input:  inode = node number
!*
!* Output:  nodes(inode)%lsq_inv_2x2 = inverse matrix for LSQ reconstruction
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
  SUBROUTINE lsq_inv_2x2(inode)
    USE constants
    USE mesh_data, ONLY : nodes

    IMPLICIT NONE
    
    INTEGER,INTENT(IN)              :: inode
    INTEGER                         :: k, isuno
    DOUBLE PRECISION                :: dx, dy, det
    DOUBLE PRECISION,DIMENSION(2,2) :: A
    
    A = 0.0d0
    
    !  Loop over the neighbor nodes.
    do k = 1,nodes(inode)%nsuno
       isuno = nodes(inode)%suno(k)
       
       dx = nodes(isuno)%x - nodes(inode)%x
       dy = nodes(isuno)%y - nodes(inode)%y
       
       A(1,1) = A(1,1) + dx*dx
       A(1,2) = A(1,2) + dx*dy
       
       A(2,1) = A(2,1) + dx*dy
       A(2,2) = A(2,2) + dy*dy
    enddo

    det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
    if (abs(det) < 1.0e-14) write(*,*) " Singular: LSQ det = ", det, " i=",inode

    ! OK, invert and store the inverse matrix:
    nodes(inode)%lsq_inv(1,1) =  A(2,2)/det
    nodes(inode)%lsq_inv(1,2) = -A(2,1)/det
    nodes(inode)%lsq_inv(2,1) = -A(1,2)/det
    nodes(inode)%lsq_inv(2,2) =  A(1,1)/det

  end subroutine lsq_inv_2x2

!!$!********************************************************************************
!!$!* Compute primative variables from conservative variables
!!$!*
!!$!* ------------------------------------------------------------------------------
!!$!*  Input:  ucons = conservative variables (d, d*vx,  d*vy, d*E)
!!$!* Output:  wprim =    primitive variables (d,    vx,   vy,   p)
!!$!* ------------------------------------------------------------------------------
!!$!* 
!!$!********************************************************************************
!!$  subroutine cons2prim(istart,iend)
!!$
!!$    USE constants, only : half,nunks,gamma,pgmin
!!$    USE mesh_data, only : ucons,wprim
!!$
!!$    IMPLICIT NONE
!!$    
!!$    INTEGER,INTENT(IN)                :: istart,iend !starting and ending node
!!$    INTEGER                           :: i,j,k
!!$    DOUBLE PRECISION                  :: KE !kinetic energy
!!$
!!$    do i=istart,iend
!!$       KE          = half*(ucons(i)%mx*ucons(i)%mx + ucons(i)%my*ucons(i)%my)/ucons(i)%d
!!$       wprim(i)%d  = ucons(i)%d
!!$       wprim(i)%vx = ucons(i)%mx/ucons(i)%d
!!$       wprim(i)%vy = ucons(i)%my/ucons(i)%d
!!$       wprim(i)%pg = max(pgmin,(gamma-1.0d0)*(ucons(i)%en - KE))
!!$    enddo
!!$
!!$  end subroutine cons2prim
!!$
!!$!********************************************************************************
!!$!* Compute conservative variables from primative variables
!!$!*
!!$!* ------------------------------------------------------------------------------
!!$!*  Input:  wprim =    primitive variables (d,    vx,   vy,   p)
!!$!* Output:  ucons = conservative variables (d, d*vx,  d*vy, d*E)
!!$!* ------------------------------------------------------------------------------
!!$!* 
!!$!********************************************************************************
!!$  subroutine prim2cons(istart,iend)
!!$
!!$    USE constants, only : half,nunks,gamma,pgmin
!!$    USE mesh_data, only : ucons,wprim
!!$
!!$    IMPLICIT NONE
!!$    
!!$    INTEGER,INTENT(IN)                :: istart,iend !starting and ending node
!!$    INTEGER                           :: i,j,k
!!$    DOUBLE PRECISION                  :: KE !kinetic energy
!!$
!!$    do i=istart,iend
!!$       KE          = half*(wprim(i)%vx*wprim(i)%vx + wprim(i)%vy*wprim(i)%vy)
!!$       ucons(i)%d  = wprim(i)%d
!!$       ucons(i)%mx = wprim(i)%d*wprim(i)%vx
!!$       ucons(i)%my = wprim(i)%d*wprim(i)%vy
!!$       ucons(i)%en = wprim(i)%pg/(gamma - 1.0d0) + wprim(i)%d*KE
!!$    enddo
!!$
!!$  end subroutine prim2cons

!********************************************************************************
!* This subroutine computes the residual norms: L1, L2, L_infty
!*
!* ------------------------------------------------------------------------------
!*  Input:  node(:)res = the residuals
!*
!* Output:  res_norm   = residual norms (L1, L2, Linf)
!* ------------------------------------------------------------------------------
!*
!* NOTE: It is not done here, but I advise you to keep the location of the
!*       maximum residual (L_inf).
!*
!********************************************************************************
  subroutine residual_norm(res_norm)

    use constants
    use mesh_data, only : nodes
    
    implicit none
    
    DOUBLE PRECISION,DIMENSION(4,3),INTENT(OUT)   :: res_norm
    !Local variables
    DOUBLE PRECISION,DIMENSION(4) :: residual
    integer :: i
    
    res_norm(:,1) =  0.0d0
    res_norm(:,2) =  0.0d0
    res_norm(:,3) = -1.0d0
    
    !--------------------------------------------------------------------------------
    do i = 1,nnode
       !--------------------------------------------------------------------------------
       residual = abs( nodes(i)%res/nodes(i)%ml )      !Divided residual
       res_norm(:,1) = res_norm(:,1)    + residual    !L1   norm
       res_norm(:,2) = res_norm(:,2)    + residual**2 !L2   norm
       res_norm(:,3) = max(res_norm(:,3), residual)   !Linf norm
       !--------------------------------------------------------------------------------
    enddo
    !--------------------------------------------------------------------------------
    
    res_norm(:,1) =      res_norm(:,1)/dble(nnode)
    res_norm(:,2) = dsqrt(res_norm(:,2)/dble(nnode))
    
  end subroutine residual_norm
  !--------------------------------------------------------------------------------

END MODULE solver
