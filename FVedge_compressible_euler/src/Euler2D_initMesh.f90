!!$ Program       : Euler2D_initMesh.f90
!!$ Author        : Andrew Kercher
!!$ Last Updated  : 02/09/2012
!!$ Description   : Module to create grid data structures
!!$ Acknowledgment: The following algorithms and diagrams were taken from 
!!$                 or heavily based off of code originally written by 
!!$                 Dr. Katate Masatsuka (info[at]cfdbooks.com).

MODULE initMesh

CONTAINS

!!$*******************************************************************************
!!$ Creates the following structures
!!$ 1. nodes: nodes of mesh
!!$     nsuel = number of surrounding elements
!!$     nsupt = number of surrounding points
!!$     ML    = lumped mass entry (a.k.a. dual volume for FV)
!!$     x,y   = point coordinates
!!$     supt  = list of surrounding points
!!$     suel  = list of surrounding elements
!!$ 2. elems: elements of mesh
!!$     nnoel = number of nodes per element
!!$     nsuel = number of surrounding elements
!!$     vol   = element volume (area in 2D)
!!$     xc,yc = element centered point coordinates
!!$     lnode = list of nodes
!!$     suel  = list of surrounding elements
!!$ 3. edges: edges of mesh
!!$     ni,nj   = nodes of edge
!!$     el1,el2 = elements connected through edge
!!$     sf      = magnitude of unit directed face vector
!!$     se      = magnitude of unit directed edgevector
!!$     nf      = unit directed face vector
!!$     ne      = unit directed edge vector
!!$ 4. bouns: boundary data for mesh
!!$     btype  = type of boundary (0:wall;
!!$     nbnode = number of boundary nodes
!!$     nbface = number of boundary faces
!!$     bnode  = list boundary nodes     
!!$     belem  =list boundary elements          
!!$     bface  = list boundary faces
!!$     nx     = x-component of outward normal                         
!!$     ny     = y-component of outward normal                         
!!$     fnx    = x-component of face outward normal
!!$     fny    = y-component of face outward normal
!!$*******************************************************************************
  SUBROUTINE initStruc(ndim)
    USE constants
    USE mesh_data
    USE solver, ONLY : lsq_inv_2x2

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim
    INTEGER            :: i,j,k
    INTEGER            :: ielem,inode,jelem,jnode,kelem,knode,iedge,jedge,kedge
    INTEGER            :: iboun,jboun,kboun,imin,imax
    INTEGER            :: suno_avg,suno_min,suno_max
    INTEGER            :: im,in
    INTEGER            :: ni,nj,nk,n1,n2,e1,e2
    INTEGER            :: nti,ntj,ntk
    DOUBLE PRECISION   :: xi,xj,xk,yi,yj,yk
    DOUBLE PRECISION   :: xc,xm,yc,ym
    LOGICAL            :: found

    !calculate bounding box
    xmin = 1.0E+16
    xmax = -xmin
    ymin = xmin
    ymax = xmax

    do i=1,nnode
       xmin = min(xmin,nodes(i)%x)
       xmax = max(xmax,nodes(i)%x)
       ymin = min(ymin,nodes(i)%y)
       ymax = max(ymax,nodes(i)%y)
    enddo

    Lx = abs(xmax - xmin)
    Ly = abs(ymax - ymin)

    !Print the data
    write(*,*) " Bounding Box:"
    write(*,'(a5,F12.8)') " xmin",xmin 
    write(*,'(a5,F12.8)') " ymin",ymin 
    write(*,'(a5,F12.8)') " xmax",xmax 
    write(*,'(a5,F12.8)') " ymax",ymax 
    write(*,*)


    !initialize to zero
    do inode=1,nnode
       nodes(inode)%nsuel = 0
       nodes(inode)%ml = 0
      nullify(nodes(inode)%suel)
    enddo

    nedge = 0

!--------------------------------------------------------------------------------
! Loop over elements and construct the fololowing data.
!
! 1. Surrounding elements: nodes(:)%nelm, nodes(:)%suel(:)
!
!    Example: Node i is surrounded by the eleemnts, 23, 101, 13, 41.
!             nodes(i)%nelm = 4
!             nodes(i)%suel(1) = 23
!             nodes(i)%suel(2) = 13
!             nodes(i)%suel(3) = 41
!             nodes(i)%suel(4) = 101
!
!        o-------o-------------o
!       /        |   .         |
!      /    23   |      41     |
!     o----------o-------------o
!      \        i \            |
!       \   101    \     13    |
!        \          \          | 
!         o----------o---------o
!
! 2. Element quantities  : elems(:)%x,elems(:)%y,elems(:)%vol
!
!  o-----------o            
!   \          |            o
!    \    (x,y)|           / \
!     \   .    |          /   \
!      \       |         /  .  \    (x,y): centroid coordinates
!       \      |        / (x,y) \     vol: volume of element
!        o-----o       o---------o



    !Elements surrounding points
    element_loop: do ielem = 1,nelem

       ni = elems(ielem)%lnode(1)
       nj = elems(ielem)%lnode(2)
       nk = elems(ielem)%lnode(3)

       xi = nodes(ni)%x
       xj = nodes(nj)%x
       xk = nodes(nk)%x

       yi = nodes(ni)%y
       yj = nodes(nj)%y
       yk = nodes(nk)%y

       !add element to list for each node
       nodes(ni)%nsuel = nodes(ni)%nsuel + 1

       call my_alloc_int_ptr(nodes(ni)%suel,nodes(ni)%nsuel)
       nodes(ni)%suel(nodes(ni)%nsuel) = ielem

       nodes(nj)%nsuel = nodes(nj)%nsuel + 1
       call my_alloc_int_ptr(nodes(nj)%suel,nodes(nj)%nsuel)
       nodes(nj)%suel(nodes(nj)%nsuel) = ielem

       nodes(nk)%nsuel = nodes(nk)%nsuel + 1
       call my_alloc_int_ptr(nodes(nk)%suel,nodes(nk)%nsuel)
       nodes(nk)%suel(nodes(nk)%nsuel) = ielem

       !calculate centroid and element area
       elems(ielem)%xc = third*(xi + xj + xk)
       elems(ielem)%yc = third*(yi + yj + yk)
       elems(ielem)%vol = tri_area(xi,xj,xk,yi,yj,yk)

       !add to lumped mass (a.k.a. dual volume)
       nodes(ni)%ml = nodes(ni)%ml + third*elems(ielem)%vol
       nodes(nj)%ml = nodes(nj)%ml + third*elems(ielem)%vol
       nodes(nk)%ml = nodes(nk)%ml + third*elems(ielem)%vol

    enddo element_loop

!--------------------------------------------------------------------------------
! Loop over elements 2
!
!  Allocate elems(:)%suel(:) : elm(:)%nsuel, elems(:)%suel(:)
!  Construct element surrounding element data: elems(:)%suel(:)
!  Order of neighbor elements [e1,e2,e3,..] are closely related to
!  the order of nodes [ni,nj,nk,..] (see below).
!
!          o------o
!          |      |                
!        nl|  e1  |nk                     nk
!    o-----o------o------o      o---------o------------o
!    |     |      |      |       .      .   .        .
!    | e2  |      |  e4  |        . e2 .     . e1  .
!    o-----o------o------o         .  .       .  .
!       ni |     .nj              ni o---------o nj   
!          | e3 .                     .   e3  .
!          |   .                        .    .
!          |  .                           . .
!          | .                             o
!          o
!

    !3 elements surround each element
    !allocate pointer for surrounding elements
    do ielem = 1,nelem
       elems(ielem)%nsuel = 3
       allocate(elems(ielem)%suel(3))
    enddo

    !elements surrounding elements
    element_loop2: do ielem = 1,nelem
       node_loop: do inode = 1,elems(ielem)%nnoel
          !get face of element
          if(inode .lt. elems(ielem)%nnoel) n2 = elems(ielem)%lnode(inode+1)
          if(inode .eq. elems(ielem)%nnoel) n2 = elems(ielem)%lnode(1)
          n1 = elems(ielem)%lnode(inode)
          !loop over elements surrounding n1
          found = .false.
          esuno: do jelem = 1,nodes(n1)%nsuel
             kelem = nodes(n1)%suel(jelem)
             do jnode = 1,elems(kelem)%nnoel
                ni = elems(kelem)%lnode(jnode)
                if(jnode .gt. 1) nj = elems(kelem)%lnode(jnode-1)
                if(jnode .eq. 1) nj = elems(kelem)%lnode(elems(kelem)%nnoel)
                if((ni .eq. n1) .and. (nj .eq. n2))then
                   found = .true.
                   im = jnode + 1
                   if(im .gt. elems(kelem)%nnoel) im = im - elems(kelem)%nnoel
                   exit 
                endif
             enddo

             if(found) exit
          enddo esuno
          
          in = inode + 2
          if(in .gt. elems(ielem)%nnoel) in = in - elems(ielem)%nnoel

          if(found) then
             elems(ielem)%suel(in) = kelem
             elems(kelem)%suel(im) = ielem
          else
             elems(ielem)%suel(in) = 0
          endif

       enddo node_loop
    enddo element_loop2

!--------------------------------------------------------------------------------
! Edge-data for node-centered (edge-based) scheme.
!
! Loop over elements 3
! Construct edge data: edges(:)%ni, nj, e1, e2.
! Edge points from node ni to node nj.
!
!      nj
!       o------------o
!     .  \         .
!    .    \   e2  .
!   .  e1  \    .
!  .        \ .        Directed area is positive: ni -> nj
! o----------o         e1: left element
!             ni       e2: right element (e2 > e1 or e2 = 0)


    !only count edge if element number is less than element number of neighbor
    do ielem=1,nelem
       if(elems(ielem)%suel(3) > ielem .or. elems(ielem)%suel(3) == 0)then
          nedge = nedge + 1
       endif
       if(elems(ielem)%suel(1) > ielem .or. elems(ielem)%suel(1) == 0)then
          nedge = nedge + 1
       endif
       if(elems(ielem)%suel(2) > ielem .or. elems(ielem)%suel(2) == 0)then
          nedge = nedge + 1
       endif
    enddo

    !allocate edge data
    allocate(edges(nedge))
    edges(:)%e1 = 0
    edges(:)%e2 = 0
    nedge = 0
    do ielem=1,nelem

       ni = elems(ielem)%lnode(1)
       nj = elems(ielem)%lnode(2)
       nk = elems(ielem)%lnode(3)

       if(elems(ielem)%suel(3) > ielem .or. elems(ielem)%suel(3) == 0)then
          nedge = nedge + 1
          edges(nedge)%ni = ni
          edges(nedge)%nj = nj
          edges(nedge)%e1 = ielem
          edges(nedge)%e2 = elems(ielem)%suel(3)
       endif
       if(elems(ielem)%suel(1) > ielem .or. elems(ielem)%suel(1) == 0)then
          nedge = nedge + 1
          edges(nedge)%ni = nj
          edges(nedge)%nj = nk
          edges(nedge)%e1 = ielem
          edges(nedge)%e2 = elems(ielem)%suel(1)
       endif
       if(elems(ielem)%suel(2) > ielem .or. elems(ielem)%suel(2) == 0)then
          nedge = nedge + 1
          edges(nedge)%ni = nk
          edges(nedge)%nj = ni
          edges(nedge)%e1 = ielem
          edges(nedge)%e2 = elems(ielem)%suel(2)
       endif

    enddo

! Loop over edges
! Construct edge vector and directed area vector.
!
! Edge vector is a simple vector pointing froom ni to nj.
! For each edge, add the scaled area vector (sav) from
! the left and right elements.
!
!              nj
!   o-----------o-----------o
!   |     sav   |  sav      |
!   |       ^   |   ^       |
!   |       |   |   |       |
!   |   c - - - m - - -c    |
!   |           |           |
!   |           |           |    m: edge midpoint
!   |           |           |    c: element centroid
!   o-----------o-----------o
!                ni
!

    do iedge = 1,nedge

       ni = edges(iedge)%ni
       nj = edges(iedge)%nj
       e1 = edges(iedge)%e1
       e2 = edges(iedge)%e2
       xm = half*(nodes(ni)%x + nodes(nj)%x)
       ym = half*(nodes(ni)%y + nodes(nj)%y)

       edges(iedge)%sav(:) = 0.0d0

       !element 1 contribution
       if(e1 > 0)then
          xc = elems(e1)%xc
          yc = elems(e1)%yc
          edges(iedge)%sav(1) = -(ym - yc)
          edges(iedge)%sav(2) = xm - xc
       endif

       !element 2 contribution
       if(e2 > 0)then
          xc = elems(e2)%xc
          yc = elems(e2)%yc
          edges(iedge)%sav(1) = edges(iedge)%sav(1) - (yc - ym)
          edges(iedge)%sav(2) = edges(iedge)%sav(2) + xc - xm
       endif

       !scaled area vector and magnitude
       edges(iedge)%sa = sqrt(edges(iedge)%sav(1)**2.0d0 + edges(iedge)%sav(2)**2.0d0)
       edges(iedge)%sav(:) = edges(iedge)%sav(:)/edges(iedge)%sa

       !scaled edge vector and magnitude
       edges(iedge)%sev(1) = nodes(nj)%x - nodes(ni)%x
       edges(iedge)%sev(2) = nodes(nj)%y - nodes(ni)%y
       edges(iedge)%se = sqrt(edges(iedge)%sev(1)**2.0d0 + edges(iedge)%sev(2)**2.0d0)
       edges(iedge)%sev(:) = edges(iedge)%sev(:)/edges(iedge)%se

    enddo

!--------------------------------------------------------------------------------
! Construct nodes surrpunding nodes:
!
!        o     o
!         \   / 
!          \ /
!     o-----*-----o
!          /|
!         / |
!        /  o        *: node in interest
!       o            o: neighbors (edge-connected)
!

    do inode = 1,nnode
       nodes(inode)%nsuno = 0
       nullify(nodes(inode)%suno)
    end do

    do iedge = 1,nedge

       ni = edges(iedge)%ni
       nj = edges(iedge)%nj
       !add nj to the neighbor list of ni
       nodes(ni)%nsuno = nodes(ni)%nsuno + 1
       call my_alloc_int_ptr(nodes(ni)%suno, nodes(ni)%nsuno)
       nodes(ni)%suno(nodes(ni)%nsuno) = nj

       !add ni to the neighbor list of nj
       nodes(nj)%nsuno = nodes(nj)%nsuno + 1
       call my_alloc_int_ptr(nodes(nj)%suno, nodes(nj)%nsuno)
       nodes(nj)%suno(nodes(nj)%nsuno) = ni

    enddo

!--------------------------------------------------------------------------------
! Boundary normal at nodes constructed by accumulating the contribution
! from each boundary face normal. This vector will be used to enforce
! the tangency condition, for example.
!
!
!        Interior domain      /
!                            o
!                  .        /
!                  .       /
! --o-------o-------------o
!           j   |  .  |   j+1
!               v  .  v
!
!        Left half added to the node j, and
!       right half added to the node j+1.
!

    ! Allocate and initialize the normal vector arrays
    do i = 1,nbcs

       allocate(bouns(i)%bnx(bouns(i)%nbnode))
       allocate(bouns(i)%bny(bouns(i)%nbnode))
       allocate(bouns(i)%bn(bouns(i)%nbnode))
       
       do j = 1,bouns(i)%nbnode
          bouns(i)%bnx(j) = 0.0d0
          bouns(i)%bny(j) = 0.0d0
          bouns(i)%bn( j) = 0.0d0
       enddo
       
    enddo

    ! Compute the outward normals
    do i = 1,nbcs
       do j = 1,bouns(i)%nbnode-1
          
          xi = nodes(bouns(i)%bnode(j  ))%x
          yi = nodes(bouns(i)%bnode(j  ))%y          
          xj = nodes(bouns(i)%bnode(j+1))%x
          yj = nodes(bouns(i)%bnode(j+1))%y

          bouns(i)%bnx(j) = bouns(i)%bnx(j) + half*( -(yi-yj) )
          bouns(i)%bny(j) = bouns(i)%bny(j) + half*(   xi-xj  )

          bouns(i)%bnx(j+1) = bouns(i)%bnx(j+1) + half*( -(yi-yj) )
          bouns(i)%bny(j+1) = bouns(i)%bny(j+1) + half*(   xi-xj  )
          
       enddo
    enddo

    ! Compute the magnitude and turn (bnx,bny) into a unit vector
    do i = 1,nbcs
       do j = 1,bouns(i)%nbnode
          bouns(i)%bn(j)  = dsqrt( bouns(i)%bnx(j)**2.0d0 + bouns(i)%bny(j)**2.0d0 )
          bouns(i)%bnx(j) = bouns(i)%bnx(j) / bouns(i)%bn(j)
          bouns(i)%bny(j) = bouns(i)%bny(j) / bouns(i)%bn(j)
       enddo
    enddo

!--------------------------------------------------------------------------------
! Boundary face data
!
!      |     Domain      |
!      |                 |
!      o--o--o--o--o--o--o  <- Boundary segment
!   j= 1  2  3  4  5  6  7
!
!   In the above case, nbnode = 7, nbface = 6
!

    do i = 1,nbcs
       bouns(i)%nbface = bouns(i)%nbnode-1

       allocate(bouns(i)%fnx(    bouns(i)%nbface   ))
       allocate(bouns(i)%fny(    bouns(i)%nbface   ))
       allocate(bouns(i)%fn(     bouns(i)%nbface   ))
       allocate(bouns(i)%belem(    bouns(i)%nbface   ))
    enddo

    do i = 1,nbcs
       do j = 1,bouns(i)%nbface
          
          xi = nodes(bouns(i)%bnode(j  ))%x
          yi = nodes(bouns(i)%bnode(j  ))%y
          xj = nodes(bouns(i)%bnode(j+1))%x
          yj = nodes(bouns(i)%bnode(j+1))%y
          
          bouns(i)%fn(j)  =  dsqrt( (xi-xj)**2.0d0 + (yi-yj)**2.0d0 )
          bouns(i)%fnx(j) = -(yi-yj)/bouns(i)%fn(j)
          bouns(i)%fny(j) =  (xi-xj)/bouns(i)%fn(j)
          
       enddo
    enddo
    
    !boundary normal at nodes: outward normal
    do i = 1,nbcs
       do j = 1,bouns(i)%nbnode-1

          xi = nodes(bouns(i)%bnode(j  ))%x
          yi = nodes(bouns(i)%bnode(j  ))%y
          xj = nodes(bouns(i)%bnode(j+1))%x
          yj = nodes(bouns(i)%bnode(j+1))%y

          bouns(i)%fn(j)  =  dsqrt( (xi-xj)**2.0d0 + (yi-yj)**2.0d0 )
          bouns(i)%fnx(j) = -(yi-yj)/bouns(i)%fn(j)
          bouns(i)%fny(j) =  (xi-xj)/bouns(i)%fn(j)
          
       enddo
    enddo


! Find element adjacent to the face: belm
!
!  NOTE: This is useful to figure out what element
!        each boundary face belongs to. Boundary flux needs
!        special weighting depending on the element.
!
!      |_________|_________|________|
!      |         |         |        | 
!      |         |         |        | 
!      |_________|_________|________|
!      |         |         |        |     <- Grid (e.g., quads)
!      |         | elemb(j)|        |
!   ---o---------o---------o--------o---  <- Boundary segment
!                 j-th face
!
! elemb(j) is the element number of the element having the j-th boundary face.
!

    do i = 1,nbcs
       do j = 1,bouns(i)%nbface

          !bface is defined by the nodes ni and nj.
          ni = bouns(i)%bnode(j  )
          nj = bouns(i)%bnode(j+1)
!!$write(*,*) ni,nj
          found = .false.
          !Find the element having the bface from the elements
          !around the node ni.
          do k = 1, nodes(ni)%nsuel
             ielem = nodes(ni)%suel(k)
             do inode = 1,elems(ielem)%nnoel
                in = inode
                im = inode+1
                if (im > elems(ielem)%nnoel) im = im - elems(ielem)%nnoel !return to 1
                nti = elems(ielem)%lnode(in)
                ntj = elems(ielem)%lnode(im)
                if (nti == ni .and. ntj == nj) then
                   found = .true.
                   exit
                endif
             enddo
             if (found) exit
          end do

          if (found) then
             bouns(i)%belem(j) = ielem
          else
             write(*,*) ' Boundary-adjacent element not found. Error...'
             write(*,*) '  boundary condition type :',bouns(i)%bname
             write(*,*) '  face =',j
             write(*,*) '  with nodes =',ni,nj             

             stop
          endif

       enddo
    enddo

!--------------------------------------------------------------------------------
! Construct least-squares matrix for node-centered schemes.
!
!        o     o
!         \   / 
!          \ /
!     o-----*-----o
!          /|
!         / |
!        /  o        *: node in interest
!       o            o: neighbors (edge-connected nghbrs)
!
! Check the number of neighbor nodes (must have at least 2 neighbors)

    write(*,*) " --- Nodes surronding nodes data:"

    suno_avg = nodes(1)%nsuno
    suno_min = nodes(1)%nsuno
    suno_max = nodes(1)%nsuno
    imin = 1
    imax = 1
    if(nodes(1)%nsuno == 2)then
       write(*,*) "--- 2 surrounding nodes for the node = ", 1
    endif

    do i = 2, nnode
       suno_avg = suno_avg + nodes(i)%nsuno
       if (nodes(i)%nsuno < suno_min) imin = i
       if (nodes(i)%nsuno > suno_max) imax = i
       suno_min = min(suno_min, nodes(i)%nsuno)
       suno_max = max(suno_max, nodes(i)%nsuno)
       if (nodes(i)%nsuno==2) then
          write(*,*) "--- 2 neighbors for the node = ", i
       endif
    end do

    write(*,*) "      suno_avg = ", suno_avg/nnode
    write(*,*) "      suno_min = ", suno_min, " at node ", imin
    write(*,*) "      suno_max = ", suno_max, " at node ", imax
    write(*,*)

    ! Now, compute the inverse of the LSQ matrix at each node.

    do i = 1,nnode
       call lsq_inv_2x2(i)
    enddo
    
    return
    
  END SUBROUTINE initStruc

!!$*******************************************************************************
!!$ Compute the area of the triangle defined by the nodes, a, b, c.
!!$
!!$              c (x3,y3)
!!$              o 
!!$             / \ 
!!$            /   \
!!$ (x1,y1) a o-----o b (x2,y2)
!!$
!!$ Nodes must be ordered counterclockwise (otherwise it gives negative area)
!!$
!!$*******************************************************************************
 function tri_area(xi,xj,xk,yi,yj,yk) result(area)
 use constants, only : half
 implicit none
 DOUBLE PRECISION,intent(in) :: xi,xj,xk,yi,yj,yk
 DOUBLE PRECISION :: area

  area = half*(xi*(yj-yk) + xj*(yk-yi) + xk*(yi-yj))

 end function tri_area

!!$*******************************************************************************
!!$ This subroutine is useful to expand or shrink integer arrays.
!!$
!!$  Array, x, will be allocated if the requested dimension is 1 (i.e., n=1)
!!$  Array, x, will be expanded to the requested dimension, n, if (n > dim(x)).
!!$  Array, x, will be shrinked to the requested dimension, n, if (n < dim(x)).
!!$
!!$*******************************************************************************
  subroutine my_alloc_int_ptr(x,n)
    implicit none
    integer, intent(in) :: n
    integer, dimension(:), pointer :: x
    integer, dimension(:), pointer :: temp
    integer :: i

    if (n <= 0) then
       write(*,*) "my_alloc_int_ptr received non-positive dimension. Stop."
       stop
    endif

    ! If not allocated, allocate and return
    if (.not.(associated(x))) then

       allocate(x(n))
       return
    endif

    ! If reallocation, create a pointer with a target of new dimension.
    allocate(temp(n))
    temp = 0

    ! (1) Expand the array dimension
    if ( n > size(x) ) then
       
       do i = 1, size(x)
          temp(i) = x(i)
       end do

       ! (2) Shrink the array dimension: the extra data, x(n+1:size(x)), discarded.
    else
       
       do i = 1, n
          temp(i) = x(i)
       end do
       
    endif
    
    ! Destroy the target of x
    deallocate(x)

    ! Re-assign the pointer
    x => temp
    
    return

  end subroutine my_alloc_int_ptr
!********************************************************************************

END MODULE initMesh
