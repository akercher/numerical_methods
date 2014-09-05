!!$ Program: Euler2D_BCs.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 02/09/2012
!!$ Description: Module includes subroutines for 2D Euler solvers.
!!$ Acknowledgment: The following algorithms and diagrams were taken from 
!!$                 or heavily based off of code originally written by 
!!$                 Dr. Katate Masatsuka (info[at]cfdbooks.com).


MODULE BCs

private                        

public :: FVedge_BCs
public :: wall_normal_mass_flux

CONTAINS


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
  subroutine FVedge_BCs
    
    use constants
    use mesh_data, only : nodes,edges,bouns,elems,ucons,wprim,gradw
    use fluxes,    only : fluxRoe2D,fluxRHLL2D,fluxHLLC2D

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
    ! Close with the boundary flux using the element-based formula that is
    ! exact for linear fluxes (See Nishikawa AIAA2010-5093 for boundary weights
    ! that ensure the linear exactness for 2D/3D elements).
    !
    !      |  Interior Domain          |
    !      |        .........          |
    !      |        .       .          |
    !      |        .       .          |
    !      o--o--o-----o---------o--o--o  <- Boundary segment
    !                  n1   |   n2
    !                       v
    !                     n12 (unit face normal vector)
    !
    ! NOTE: We visit each boundary face, defined by the nodes n1 and n2,
    !       and compute the flux across the boundary face: left half for node1,
    !       and the right half for node2. In the above figure, the dots indicate
    !       the control volume around the node n1. Clearly, the flux across the
    !       left half of the face contributes to the node n1. Similarly for n2.
    !
    !
    !--------------------------------------------------------------------------------
    bc_loop : do i = 1,nbcs
       !--------------------------------------------------------------------------------

       !------------------------------------------------
       !  BC: Upwind flux via freestream values
       !
       !      NOTE: If the final solution at the boundary node is far from
       !            the freestream values, then the domain is probably is not large enough.
       
       bc : if (trim(bouns(i)%bname) == "freestream") then
          
          do j = 1, bouns(i)%nbface
             
             nfl     = bouns(i)%bnode(j)   !left node
             nfr     = bouns(i)%bnode(j+1) !right node
             bfn(1)  = bouns(i)%fnx(j)     !x-component of the unit face normal vector
             bfn(2)  = bouns(i)%fny(j)     !y-component of the unit face normal vector
             bfn_mag = half*bouns(i)%fn(j) !half length of the boundary face, j.
             
             !left node
             wl(1) = wprim(nfl)%d 
             wl(2) = wprim(nfl)%vx 
             wl(3) = wprim(nfl)%vy 
             wl(4) = wprim(nfl)%pg
             
             wr(1) = dinf
             wr(2) = vxinf
             wr(3) = vyinf
             wr(4) = pginf

             if     (trim(flux_type)=="roe") then
          
                call fluxRoe2D(bfn,wl,wr,bfl,nwsd)
                nodes(nfl)%nwsd = nodes(nfl)%nwsd + nwsd*bfn_mag
          
             elseif (trim(flux_type)=="rhll") then
                
                call fluxRHLL2D(bfn,wl,wr,bfl,nwsd)
                nodes(nfl)%nwsd = nodes(nfl)%nwsd + nwsd*bfn_mag

             elseif (trim(flux_type)=="hllc") then
                
                call fluxHLLC2D(bfn,wl,wr,bfr,nwsd)
                nodes(nfr)%nwsd = nodes(nfr)%nwsd + nwsd*bfn_mag
                
             else
          
                write(*,*) " Invalid input for flux_type = ", trim(flux_type)
                write(*,*) " Choose roe, rhll or hllc, and try again."
                write(*,*) " ... Stop."
                stop
                
             endif
             
             !right node
             wl(1) = wprim(nfr)%d 
             wl(2) = wprim(nfr)%vx 
             wl(3) = wprim(nfr)%vy 
             wl(4) = wprim(nfr)%pg
             
             wr(1) = dinf
             wr(2) = vxinf
             wr(3) = vyinf
             wr(4) = pginf
             
             if     (trim(flux_type)=="roe") then
          
                call fluxRoe2D(bfn,wl,wr,bfr,nwsd)
                nodes(nfr)%nwsd = nodes(nfr)%nwsd + nwsd*bfn_mag

             elseif (trim(flux_type)=="rhll") then
                
                call fluxRHLL2D(bfn,wl,wr,bfr,nwsd)
                nodes(nfr)%nwsd = nodes(nfr)%nwsd + nwsd*bfn_mag

             elseif (trim(flux_type)=="hllc") then
                
                call fluxHLLC2D(bfn,wl,wr,bfr,nwsd)
                nodes(nfr)%nwsd = nodes(nfr)%nwsd + nwsd*bfn_mag
                
             else
          
                write(*,*) " Invalid input for flux_type = ", trim(flux_type)
                write(*,*) " Choose roe, rhll or hllc, and try again."
                write(*,*) " ... Stop."
                stop
                
             endif
             
             !3. Add contributions to the two nodes (See Nishikawa AIAA2010-5093)
             belem = bouns(i)%belem(j)

             nodes(nfl)%res(:) = nodes(nfl)%res(:) + sixth*bfn_mag*(5.0d0*bfl(:) + bfr(:))
             nodes(nfr)%res(:) = nodes(nfr)%res(:) + sixth*bfn_mag*(5.0d0*bfr(:) + bfl(:))
!!$             nodes(nfl)%res(:) = nodes(nfl)%res(:) + bfn_mag*(5.0d0*bfl(:) + bfr(:))/6.0d0
!!$             nodes(nfr)%res(:) = nodes(nfr)%res(:) + bfn_mag*(5.0d0*bfr(:) + bfl(:))/6.0d0

          end do

       !------------------------------------------------
       !  BC: Solid body and Supersonic outflow
       !
       !      NOTE: Basically, simply compute the physical flux, which
       !            can be done by calling the Roe flux with wR = wL.
       !            It is equivalent to the interior-extrapolation condition.
       !
       !      NOTE: Tangency condition for solid body will be applied later.
       elseif (trim(bouns(i)%bname) == "slip_wall" .or.    &
            trim(bouns(i)%bname) == "outflow_supersonic") then
          
          bnodes_slip_wall : do j = 1, bouns(i)%nbface
             
             nfl     = bouns(i)%bnode(j)   !left node
             nfr     = bouns(i)%bnode(j+1) !right node
             bfn(1)  = bouns(i)%fnx(j)     !x-component of the unit face normal vector
             bfn(2)  = bouns(i)%fny(j)     !y-component of the unit face normal vector
             bfn_mag = half*bouns(i)%fn(j) !half length of the boundary face, j.
             
             !left node
             wl(1) = wprim(nfl)%d 
             wl(2) = wprim(nfl)%vx 
             wl(3) = wprim(nfl)%vy 
             wl(4) = wprim(nfl)%pg
             
             wr(:) = wl(:)
             
             if     (trim(flux_type)=="roe") then
          
                call fluxRoe2D(bfn,wl,wr,bfl,nwsd)
                nodes(nfl)%nwsd = nodes(nfl)%nwsd + nwsd*bfn_mag
          
             elseif (trim(flux_type)=="rhll") then
                
                call fluxRHLL2D(bfn,wl,wr,bfl,nwsd)
                nodes(nfl)%nwsd = nodes(nfl)%nwsd + nwsd*bfn_mag

             elseif (trim(flux_type)=="hllc") then
                
                call fluxHLLC2D(bfn,wl,wr,bfr,nwsd)
                nodes(nfr)%nwsd = nodes(nfr)%nwsd + nwsd*bfn_mag
                
             else
          
                write(*,*) " Invalid input for flux_type = ", trim(flux_type)
                write(*,*) " Choose roe, rhll or hllc, and try again."
                write(*,*) " ... Stop."
                stop
                
             endif
             
             !right node
             wl(1) = wprim(nfr)%d 
             wl(2) = wprim(nfr)%vx 
             wl(3) = wprim(nfr)%vy 
             wl(4) = wprim(nfr)%pg
             
             wr(:) = wl(:)
             
             if     (trim(flux_type)=="roe") then
          
                call fluxRoe2D(bfn,wl,wr,bfr,nwsd)
                nodes(nfr)%nwsd = nodes(nfr)%nwsd + nwsd*bfn_mag
          
             elseif (trim(flux_type)=="rhll") then
                
                call fluxRHLL2D(bfn,wl,wr,bfr,nwsd)
                nodes(nfr)%nwsd = nodes(nfr)%nwsd + nwsd*bfn_mag

             elseif (trim(flux_type)=="hllc") then
                
                call fluxHLLC2D(bfn,wl,wr,bfr,nwsd)
                nodes(nfr)%nwsd = nodes(nfr)%nwsd + nwsd*bfn_mag
                
             else
          
                write(*,*) " Invalid input for flux_type = ", trim(flux_type)
                write(*,*) " Choose roe, rhll or hllc, and try again."
                write(*,*) " ... Stop."
                stop
                
             endif
             
             !3. Add contributions to the two nodes (See Nishikawa AIAA2010-5093)
             belem = bouns(i)%belem(j)
             
             nodes(nfl)%res(:) = nodes(nfl)%res(:) + sixth*bfn_mag*(5.0d0*bfl(:) + bfr(:))
             nodes(nfr)%res(:) = nodes(nfr)%res(:) + sixth*bfn_mag*(5.0d0*bfr(:) + bfl(:))
!!$             nodes(nfl)%res(:) = nodes(nfl)%res(:) + bfn_mag*(5.0d0*bfl(:) + bfr(:))/6.0d0
!!$             nodes(nfr)%res(:) = nodes(nfr)%res(:) + bfn_mag*(5.0d0*bfr(:) + bfl(:))/6.0d0
             
          end do bnodes_slip_wall

       !------------------------------------------------
       !  BC: Subsonic Outflow - Fixed Back Pressure
       !
       !      NOTE: Fix the pressure as freestream pressure
       !            on the right side of the face (outside the domain).
       !            Assumption is that the outflow boundary is far from the body.          
       elseif (trim(bouns(i)%bname) == "outflow_back_pressure") then
          
          bnodes_outflow : do j = 1, bouns(i)%nbface
             
             nfl     = bouns(i)%bnode(j)   !left node
             nfr     = bouns(i)%bnode(j+1) !right node
             bfn(1)  = bouns(i)%fnx(j)     !x-component of the unit face normal vector
             bfn(2)  = bouns(i)%fny(j)     !y-component of the unit face normal vector
             bfn_mag = half*bouns(i)%fn(j) !half length of the boundary face, j.
             
             !left node
             wl(1) = wprim(nfl)%d 
             wl(2) = wprim(nfl)%vx 
             wl(3) = wprim(nfl)%vy 
             wl(4) = wprim(nfl)%pg
             
             wr(:) = wl(:)
             wr(4) = pginf
             
             if     (trim(flux_type)=="roe") then
          
                call fluxRoe2D(bfn,wl,wr,bfl,nwsd)
                nodes(nfl)%nwsd = nodes(nfl)%nwsd + nwsd*bfn_mag
          
             elseif (trim(flux_type)=="rhll") then
                
                call fluxRHLL2D(bfn,wl,wr,bfl,nwsd)
                nodes(nfl)%nwsd = nodes(nfl)%nwsd + nwsd*bfn_mag

             elseif (trim(flux_type)=="hllc") then
                
                call fluxHLLC2D(bfn,wl,wr,bfr,nwsd)
                nodes(nfr)%nwsd = nodes(nfr)%nwsd + nwsd*bfn_mag
                
             else
          
                write(*,*) " Invalid input for flux_type = ", trim(flux_type)
                write(*,*) " Choose roe, rhll or hllc, and try again."
                write(*,*) " ... Stop."
                stop
                
             endif
             
             !right node
             wl(1) = wprim(nfr)%d 
             wl(2) = wprim(nfr)%vx 
             wl(3) = wprim(nfr)%vy 
             wl(4) = wprim(nfr)%pg
             
             wr(:) = wl(:)
             wr(4) = pginf
                          
             if     (trim(flux_type)=="roe") then
          
                call fluxRoe2D(bfn,wl,wr,bfr,nwsd)
                nodes(nfr)%nwsd = nodes(nfr)%nwsd + nwsd*bfn_mag
          
             elseif (trim(flux_type)=="rhll") then
                
                call fluxRHLL2D(bfn,wl,wr,bfr,nwsd)
                nodes(nfr)%nwsd = nodes(nfr)%nwsd + nwsd*bfn_mag
                
             elseif (trim(flux_type)=="hllc") then
                
                call fluxHLLC2D(bfn,wl,wr,bfr,nwsd)
                nodes(nfr)%nwsd = nodes(nfr)%nwsd + nwsd*bfn_mag

             else
          
                write(*,*) " Invalid input for flux_type = ", trim(flux_type)
                write(*,*) " Choose roe, rhll or hllc, and try again."
                write(*,*) " ... Stop."
                stop
                
             endif

             !3. Add contributions to the two nodes (See Nishikawa AIAA2010-5093)
             belem = bouns(i)%belem(j)
             
             nodes(nfl)%res(:) = nodes(nfl)%res(:) + sixth*bfn_mag*(5.0d0*bfl(:) + bfr(:))
             nodes(nfr)%res(:) = nodes(nfr)%res(:) + sixth*bfn_mag*(5.0d0*bfr(:) + bfl(:))
!!$             nodes(nfl)%res(:) = nodes(nfl)%res(:) + bfn_mag*(5.0d0*bfl(:) + bfr(:))/6.0d0
!!$             nodes(nfr)%res(:) = nodes(nfr)%res(:) + bfn_mag*(5.0d0*bfr(:) + bfl(:))/6.0d0
             
          end do bnodes_outflow
          !------------------------------------------------

       !------------------------------------------------
       !  BC: Periodic
       !
       !      NOTE: Add residual from from p2 to residual of p1, 
       !            then set residual of p2 = residual of p1
       !
       !               p1 o-------------o p2
       !   

       elseif (trim(bouns(i)%bname) == "periodic") then
          
          bnodes_periodic : do j = 1, bouns(i)%nbnode
             
             ni      = bouns(i)%bnode(j)   !boundary node
             nj     = bouns(i)%bnode(bouns(i)%nbnode + j) !periodic boundary node
                          

             nodes(ni)%res(:) = nodes(ni)%res(:) + nodes(nj)%res(:)
             nodes(ni)%nwsd   = nodes(ni)%nwsd + nodes(nj)%nwsd

             nodes(nj)%res(:) = nodes(ni)%res(:)
             nodes(nj)%nwsd   = nodes(ni)%nwsd

          end do bnodes_periodic
          !------------------------------------------------
          
       endif bc
       
       !--------------------------------------------------------------------------------
    end do bc_loop
    !--------------------------------------------------------------------------------
    
    !------------------------------------------------
    !  BC: Solid body - Slip condition (Tangency condition).
    !
    !  NOTE: It is good to enforce this condition after everything else has been done.
    
    bc_loop2 : do i = 1,nbcs
       
       only_slip_wall : if (trim(bouns(i)%bname) == "slip_wall") then
          
          bnodes_slip_wall2 : do j = 1, bouns(i)%nbnode
             
             ni = bouns(i)%bnode(j)
             nbd(1) = bouns(i)%bnx(j)
             nbd(2) = bouns(i)%bny(j)
             
             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! THIS IS A SPECIAL TREATMENT FOR SHOCK DIFFRACTION PROBLEM.
             ! Same as in the subroutine "eliminate_normal_mass_flux" below.
             if(trim(prob_type) == 'shock_diff')then
                if (i==2 .and. j==1) then

                   nodes(ni)%res(3) = 0.0d0  !no updates to y-momentum.

                   cycle bnodes_slip_wall2   !Go to the next.
                endif
             endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
             !    Subtract the normal component of the mass flow for tangency, so that
             !    normal mass flow through the boundary will not be created at nodes
             !    after the solution update.
             
             mnorm = nodes(ni)%res(2)*nbd(1) + nodes(ni)%res(3)*nbd(2)
             nodes(ni)%res(2) = nodes(ni)%res(2) - mnorm*nbd(1)
             nodes(ni)%res(3) = nodes(ni)%res(3) - mnorm*nbd(2)
             
          end do bnodes_slip_wall2
          
       endif only_slip_wall
       
    end do bc_loop2
        
  end subroutine FVedge_BCs
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
  subroutine wall_normal_mass_flux
    
    use constants, only : nunks,nbcs,prob_type
    use mesh_data, only : nodes,bouns,ucons,wprim
    use state_variables

    IMPLICIT NONE

    INTEGER                       :: i,j,inode
    DOUBLE PRECISION              :: nmf     !normal mass flux
    DOUBLE PRECISION,DIMENSION(2) :: nbd     !boundary normal

    bc_loop : do i = 1, nbcs

       only_slip_wall : if(trim(bouns(i)%bname) == "slip_wall") then

          write(*,*) " Eliminating the normal momentum on slip wall boundary ", i

          bnode_slip_wall : do j = 1,bouns(i)%nbnode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS IS A SPECIAL TREATMENT FOR SHOCK DIFFRACTION PROBLEM.
!
! NOTE: This is a corner point between the inflow boundary and
!       the lower-left wall. Enforce zero y-momentum, which is
!       not ensured by the standard BCs.
!       This special treatment is necessary because the domain
!       is rectangular (the left boundary is a straight ine) and
!       the midpoint node on the left boundary is actually a corner.
!
!       Our computational domain:
!
!                 ---------------
!          Inflow |             |
!                 |             |  o: Corner node
!          .......o             |
!            Wall |             |  This node is a corner.
!                 |             |
!                 ---------------
!
!       This is to simulate the actual domain shown below:
!      
!         -----------------------
! Inflow  |                     |
!         |                     |  o: Corner node
!         --------o             |
!            Wall |             |
!                 |             |
!                 ---------------
!      In effect, we're simulating this flow by a simplified
!      rectangular domain (easier to generate the grid).
!      So, an appropriate slip BC at the corner node needs to be applied,
!      which is "zero y-momentum", and that's all.
!
             if(trim(prob_type) == 'shock_diff')then
                if (i==2 .and. j==1) then
                   inode = bouns(i)%bnode(j)
                   ucons(inode)%my = 0.0d0     !set my to zero
                   call cons2prim(inode,inode) !update primitive variables
                   cycle bnode_slip_wall       !go to the next.
                endif
             endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             inode = bouns(i)%bnode(j)
             nbd(1) = bouns(i)%bnx(j)
             nbd(2) = bouns(i)%bny(j)
             
             !inner product of momentum and normal
             nmf = ucons(inode)%mx*nbd(1) + ucons(inode)%my*nbd(2)

             !eliminate normal mass flux from boundary
             ucons(inode)%mx = ucons(inode)%mx - nmf*nbd(1)
             ucons(inode)%my = ucons(inode)%my - nmf*nbd(2)
             
             call cons2prim(inode,inode) !Update primitive variables

          end do bnode_slip_wall

          write(*,*) " Finished eliminating the normal momentum on slip wall boundary ", i
          write(*,*)
          
       endif only_slip_wall
       
    end do bc_loop

  end subroutine wall_normal_mass_flux
!--------------------------------------------------------------------------------

END MODULE BCs
