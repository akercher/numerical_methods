!!$ Program       : FV_mesh_check.f90
!!$ Author        : Andrew Kercher
!!$ Last Updated  : 02/09/2012
!!$ Description   : Module to check mesh data for finite volume edge-based
!!$                 Euler solver.
!!$ Acknowledgment: The following algorithms and diagrams were taken from 
!!$                 or heavily based off of code originally written by 
!!$                 Dr. Katate Masatsuka (info[at]cfdbooks.com).

MODULE mesh_check

CONTAINS

!********************************************************************************
!* Check the grid data.
!*
!* 1. Directed area must sum up to zero around every node.
!* 2. Directed area must sum up to zero over the entire grid.
!* 3. Global sum of the boundary normal vectors must vanish.
!* 4. Global sum of the boundary face normal vectors must vanish.
!* 5. Check element volumes which must be positive.
!* 6. Check dual volumes which must be positive.
!* 7. Global sum of the dual volumes must be equal to the sum of element volumes.
!* 8. Linear LSQ gradients must be exact for linear functions.
!*
!* Add more tests you can think of.
!*
!********************************************************************************
  SUBROUTINE mesh_tests

    USE constants
    USE mesh_data, ONLY : nnode,nelem,nedge,nodes,elems,nedge,edges,nbcs,bouns,wprim,gradw
    USE solver, ONLY : lsq_grads_2x2

    IMPLICIT NONE
    !Local variables
    INTEGER                              :: i,j,k
    INTEGER                              :: ni,nj,ierr
    DOUBLE PRECISION                     :: sum_volc, sum_vol
    DOUBLE PRECISION                     :: c0_linear, cx_linear, cy_linear
    DOUBLE PRECISION, DIMENSION(nnode,2) :: sum_sav_i
    DOUBLE PRECISION, DIMENSION(2)       :: sum_sav, sum_bn
    DOUBLE PRECISION, DIMENSION(2)       :: sum_bfn
    
    write(*,*) "Checking grid data...."

!--------------------------------------------------------------------------------
! Directed area sum check
!--------------------------------------------------------------------------------

    ! Compute the sum of the directed area for each node.
    sum_sav_i = 0.0d0
    do i = 1,nedge
       ni = edges(i)%ni
       nj = edges(i)%nj
       sum_sav_i(ni,:) = sum_sav_i(ni,:) + edges(i)%sav(:)*edges(i)%sa
       sum_sav_i(nj,:) = sum_sav_i(nj,:) - edges(i)%sav(:)*edges(i)%sa
    end do

    ! Compute also the sum of the boundary normal vector (at nodes).
    sum_bn = 0.0d0
    do i = 1, nbcs
       do j = 1, bouns(i)%nbnode
          k = bouns(i)%bnode(j)
          sum_sav_i(k,1) = sum_sav_i(k,1) + bouns(i)%bnx(j)*bouns(i)%bn(j)
          sum_sav_i(k,2) = sum_sav_i(k,2) + bouns(i)%bny(j)*bouns(i)%bn(j)
          sum_bn(1) = sum_bn(1) + bouns(i)%bnx(j)*bouns(i)%bn(j)
          sum_bn(2) = sum_bn(2) + bouns(i)%bny(j)*bouns(i)%bn(j)
       enddo
    enddo

    ! Global sum of boundary normal vectors must vanish.
    if (sum_bn(1) > 1.0e-12 .and. sum_bn(2) > 1.0e-12) then
       write(*,*) "--- Global sum of the boundary normal vector:"
       write(*,'(a19,es10.3)') "    sum of bn_x = ", sum_bn(1)
       write(*,'(a19,es10.3)') "    sum of bn_y = ", sum_bn(2)
       write(*,*) "Error: boundary normal vectors do not sum to zero..."
       stop
    endif

    ! Sum of the directed area vectors must vanish at every node.
    do i = 1,nnode
       if (abs(sum_sav_i(i,1))>1.0e-12 .or. abs(sum_sav_i(i,2))>1.0e-12) then
          write(*,'(a11,i5,a7,2es10.3,a9,2es10.3)') &
               " --- node=", i, " (x,y)=", nodes(i)%x, nodes(i)%y, " sum_sav=",sum_sav_i(i,:)
       endif
    enddo
    
    write(*,*) "--- Max sum of directed area vector around a node:"
    write(*,*) "  max(sum_sav_i_x) = ", maxval(sum_sav_i(:,1))
    write(*,*) "  max(sum_sav_i_y) = ", maxval(sum_sav_i(:,2))
    
    if (maxval(abs(sum_sav_i(:,1)))>1.0e-12 .or. &
         maxval(abs(sum_sav_i(:,2)))>1.0e-12) then
       write(*,*) "--- Max sum of directed area vector around a node:"
       write(*,*) "  max(sum_sav_i_x) = ", maxval(sum_sav_i(:,1))
       write(*,*) "  max(sum_sav_i_y) = ", maxval(sum_sav_i(:,2))
       write(*,*) "Error: directed area vectors do not sum to zero..."
       stop
    endif
    
    ! Of course, the global sum of the directed area vector sum must vanish.
    sum_sav = 0.0d0
    do i = 1, nnode
       sum_sav = sum_sav + sum_sav_i(i,:)
    enddo
    
    write(*,*) "--- Global sum of the directed area vector:"
    write(*,'(a19,es10.3)') "    sum of sav_x = ", sum_sav(1)
    write(*,'(a19,es10.3)') "    sum of sav_y = ", sum_sav(2)
    
    if (sum_sav(1) > 1.0e-12 .and. sum_sav(2) > 1.0e-12) then
       write(*,*) "Error: directed area vectors do not sum globally to zero..."
       write(*,*) "--- Global sum of the directed area vector:"
       write(*,'(a19,es10.3)') "    sum of dav_x = ", sum_sav(1)
       write(*,'(a19,es10.3)') "    sum of dav_y = ", sum_sav(2)
       stop
    endif
    
!--------------------------------------------------------------------------------
! Global sum check for boundary face vector
!--------------------------------------------------------------------------------
  sum_bfn(:) = 0.0d0
  do i = 1, nbcs
   do j = 1, bouns(i)%nbface
     sum_bfn(1) =  sum_bfn(1) + bouns(i)%fnx(j)*bouns(i)%fn(j)
     sum_bfn(2) =  sum_bfn(2) + bouns(i)%fny(j)*bouns(i)%fn(j)
   end do
  end do

   write(*,*) "--- Global sum of the boundary face vector:"
   write(*,'(a19,es10.3)') "    sum of bfn_x = ", sum_bfn(1)
   write(*,'(a19,es10.3)') "    sum of bfn_y = ", sum_bfn(2)

  if (sum_bfn(1) > 1.0e-12 .and. sum_bfn(2) > 1.0e-12) then
   write(*,*) "Error: boundary face normals do not sum globally to zero..."
   write(*,*) "--- Global sum of the boundary face normal vector:"
   write(*,'(a19,es10.3)') "    sum of bfn_x = ", sum_bfn(1)
   write(*,'(a19,es10.3)') "    sum of bfn_y = ", sum_bfn(2)
   stop
  endif

!--------------------------------------------------------------------------------
! Volume check
!--------------------------------------------------------------------------------
! (1)Check the element volume: make sure there are no zero or negative volumes
  ierr = 0
  sum_volc = 0.0d0
  do i = 1,nelem

     sum_volc = sum_volc + elems(i)%vol
     
     if (elems(i)%vol < 0.0d0) then
        write(*,*) "Negative volc=",elems(i)%vol, " elem=",i, " stop..."
        ierr = ierr + 1
     endif

     if (dabs(elems(i)%vol) < 1.0e-14) then
        write(*,*) "Vanishing volc=",elems(i)%vol, " elem=",i, " stop..."
        ierr = ierr + 1
     endif
  end do

!--------------------------------------------------------------------------------
! (2)Check the mass matraix (a.k.a. dual volume) around a node
  ierr = 0
  sum_vol = 0.0d0
  do i = 1,nnode
     
     sum_vol = sum_vol + nodes(i)%ml

     if (nodes(i)%ml < 0.0d0) then
        write(*,*) "Negative vol=",nodes(i)%ml, " node=",i, " stop..."
        ierr = ierr + 1
     endif

     if (abs(nodes(i)%ml) < tol) then
        write(*,*) "Vanishing vol=",nodes(i)%ml, " node=",i, " stop..."
        ierr = ierr + 1
     endif

  end do

  if (ierr > 0) stop

  !sum of global volume must equal sum of mass matrix
  if (abs(sum_vol-sum_volc) > 1.0e-11) then
     write(*,*) "--- Global sum of volume: must be the same"
     write(*,'(a19,es10.3)') "    sum of volc = ", sum_volc
     write(*,'(a19,es10.3)') "    sum of vol  = ", sum_vol
     write(*,'(a22,es10.3)') " sum_vol-sum_volc  = ", sum_vol-sum_volc
     write(*,*) "Error: sum of dual volumes and cell volumes do not match..."
     stop
  endif
  
!--------------------------------------------------------------------------------
! Check the least-squares matrix for node-centered scheme.
! Check if the LSQ gradients are exact for a linear function.
!--------------------------------------------------------------------------------
  write(*,*) "Checking the least-squares matrix(nc,linear)..."

  ! Store a linear function
  c0_linear = 5.0d0
  cx_linear = 3.70d0
  cy_linear = 2.30d0

  !allocate primative variables and gradients
  allocate(wprim(nnode),gradw(nnode))

  do i = 1,nnode
     wprim(i)%d  = cx_linear*nodes(i)%x + cy_linear*nodes(i)%y + c0_linear
     wprim(i)%vx = cx_linear*nodes(i)%x + cy_linear*nodes(i)%y + c0_linear
     wprim(i)%vy = cx_linear*nodes(i)%x + cy_linear*nodes(i)%y + c0_linear
     wprim(i)%pg = cx_linear*nodes(i)%x + cy_linear*nodes(i)%y + c0_linear
  end do
    
!-----------------------------------------------------------------
! Check the 2x2 gradients for node-centered schemes.
!-----------------------------------------------------------------
! Compute gradients

  do i = 1, nnode
     call lsq_grads_2x2(i)
  end do

  ! Check the gradients
  do i = 1,nnode

     !  Check dw/dx
     if (dabs( maxval(gradw(i)%dx(:)) - cx_linear) > 1.0e-11) then
        write(*,*) " i = ", i, "  nsuno=",nodes(i)%nsuno
        write(*,*) " Max error = ", maxval(gradw(i)%dx(:)) - cx_linear
        write(*,'(a9,9es10.2)') " gradw_x=", gradw(i)%dx(:)
        stop
     endif

     !  Check dw/dy
     if (dabs( maxval(gradw(i)%dy(:)) - cy_linear) > 1.0e-11) then
        write(*,*) " i = ", i, "  nsuno=",nodes(i)%nsuno
        write(*,*) " Max error = ", maxval(gradw(i)%dy(:)) - cy_linear
        write(*,'(a9,9es10.2)') " gradw_y=", gradw(i)%dy(:)
        stop
     endif

  end do

  !nullify primative variables and their gradients
  nullify(wprim,gradw)

  write(*,*)
  write(*,*) "No problems detected for grid data.  Proceeding with calculation."  

  END SUBROUTINE mesh_tests


END MODULE mesh_check
