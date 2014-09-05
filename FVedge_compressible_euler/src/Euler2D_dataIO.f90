!!$ Program:      Euler2D_dataIO.f90
!!$ Author:       Andrew Kercher
!!$ Last Updated: 02/09/2012
!!$ Description:  Module to handle input data and write output data.

MODULE dataIO

  INTERFACE dataInput
     MODULE PROCEDURE dataInput_easy,dataInput_mpfem
  END INTERFACE

CONTAINS

!********************************************************************************
!* read input file
  SUBROUTINE dataInput_easy(data_mesh,data_bcs)
    USE constants
    USE data_types
    USE mesh_data

    IMPLICIT NONE

    INTEGER :: i,j,k,os
    INTEGER :: bufferI          !buffer for integer values
    CHARACTER (LEN=50),INTENT(IN) :: data_mesh,data_bcs
    CHARACTER (LEN=30) :: datasetName,fldName,meshName,endflds
    CHARACTER (LEN=50):: bufferC ! Buffer for information that is disregarded
    
    ! Open file containing input data
    open(unit=1001, file = data_mesh, status="unknown", iostat=os)

    ! READ: Get the size of the grid.
    read(1001,*) nnode,nelem,bufferI

    !  Allocate node and element arrays.
    allocate(nodes(nnode))
    allocate(elems(nelem))

    ! READ: Read coordinates of nodes
    do i = 1, nnode
       read(1001,*) nodes(i)%x, nodes(i)%y
    enddo

    ! Read element-connectivity information

    ! Triangles: assumed that the vertices are ordered counterclockwise
    !
    !         v3
    !         /\
    !        /  \
    !       /    \
    !      /      \
    !     /        \
    !    /__________\
    !   v1           v2

    do i = 1,nelem
       elems(i)%nnoel = 3
       allocate(elems(i)%lnode(3))
       read(1001,*) elems(i)%lnode(1),elems(i)%lnode(2), elems(i)%lnode(3)
    end do
    

    !  Write out the grid data.
    
    write(*,*) " Total numbers:"
    write(*,*) "      nodes = ", nnode
    write(*,*) "   elements = ", nelem
    write(*,*)
    
    ! Read the boundary grid data
    
    ! READ: Number of boundary condition types
    read(1001,*) nbcs
    allocate(bouns(nbcs))

    ! 2. Read the boundary condition data file
    write(*,*) "Reading the boundary condition file....", data_bcs

    ! Open the input file.
    open(unit=1002, file = data_bcs, status="unknown", iostat=os)

    read(1002,*) 
    ! READ: Read the boundary condition type
    do i = 1,nbcs
       print*,i
       read(1002,*) bufferI, bouns(i)%bname
    end do

    ! READ: Number of Boundary nodes (including the starting one at the end if
    ! it is closed such as an airfoil.)
    do i = 1,nbcs
       read(1001,*) bouns(i)%nbnode
       allocate(bouns(i)%bnode(bouns(i)%nbnode))
    enddo
    
    ! READ: Read boundary nodes
    do i = 1,nbcs
       do j = 1,bouns(i)%nbnode
          read(1001,*) bouns(i)%bnode(j)
       end do
       !reduce number of boundary nodes for periodic boundaries
       if(trim(bouns(i)%bname) == 'periodic')then
          bouns(i)%nbnode = bouns(i)%nbnode/2
       endif
    end do

    !  Print the boundary grid data.
    write(*,*) " Boundary nodes:"
    write(*,*) "    segments = ", nbcs
    do i = 1, nbcs
       write(*,'(a9,i3,2(a11,i5))') " boundary", i, "  bnodes = ", bouns(i)%nbnode, &
            "  bfaces = ", bouns(i)%nbnode-1
    end do
    write(*,*)
    
    close(1001)


    !  Print the data
    write(*,*) " Boundary conditions:"
    do i = 1, nbcs
       write(*,'(a9,i3,a12,a35)') " boundary", i, "  bc_type = ", trim(bouns(i)%bname)
    end do

    write(*,*)
    
    close(1002)

  END SUBROUTINE dataInput_easy
!********************************************************************************

!********************************************************************************
!* read input file
  SUBROUTINE dataInput_mpfem(ndim,fname)
    USE constants
    USE data_types
    USE mesh_data

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim
    INTEGER :: i,j,k
    INTEGER :: bufferI          !buffer for integer values
    CHARACTER (LEN=50),INTENT(IN) :: fname
    CHARACTER (LEN=30) :: datasetName,fldName,meshName,endflds
    CHARACTER (LEN=50):: bufferC ! Buffer for information that is disregarded
    
    ! Open file containing input data
    OPEN(1001,FILE = fname)

    READ(1001,*) bufferI
    READ(1001,*) bufferC        
    READ(1001,*) bufferC        
    READ(1001,*) bufferC        
    READ(1001,*) bufferC        !ndimn ntype
    READ(1001,*) bufferI
    READ(1001,*) bufferC        !nelem npoin nboun time
    READ(1001,*) bufferI
    READ(1001,*) bufferC        

    !read in connectivity array
    DO i=1,nelem
       elems(i)%nnoel = 3
       ALLOCATE(elems(i)%lnode(3))
       READ(1001,*) bufferI,elems(i)%lnode(1),elems(i)%lnode(2),elems(i)%lnode(3)
    ENDDO

    READ(1001,*) bufferC
    !read in point coordinates
    DO i=1,nnode
       READ(1001,*) bufferI,nodes(i)%x,nodes(i)%y
    ENDDO

    !read boundary points
    !boundary type, 0:wall, 4:farfield
    !read in number of types of BCs
    READ(1001,*) bufferC
    READ(1001,*) nbcs
    ALLOCATE(bouns(nbcs))

    !read number of boundary points for each condition 
    !add extra boundary node equal to first boundary node if 
    !boundary is closed
    DO i = 1,nbcs
       READ(1001,*) bouns(i)%nbnode, bouns(i)%btype
       !if all the same condition then boundary is closed
       IF(nbcs == 1)THEN
          bouns(i)%nbnode = bouns(i)%nbnode + 1          
       !if airfoil add extra node to end that is equal to first node
       ELSEIF(bouns(i)%btype .EQ. 0)THEN
          bouns(i)%nbnode = bouns(i)%nbnode + 1
       ENDIF
       ALLOCATE(bouns(i)%bnode(bouns(i)%nbnode))
    ENDDO

    DO i = 1,nbcs
       DO j = 1,bouns(i)%nbnode-1
          READ(1001,*) bouns(i)%bnode(j),bufferI
       ENDDO
       !if closed boundary add extra node to end that is equal to first node
       IF(nbcs .EQ. 1)THEN
          bouns(i)%bnode(bouns(i)%nbnode) = bouns(i)%bnode(1)
       ELSEIF(bouns(i)%btype .EQ. 0)THEN
          bouns(i)%bnode(bouns(i)%nbnode) = bouns(i)%bnode(1)
       ELSE
          !read in last node
          READ(1001,*) bouns(i)%bnode(j),bufferI
       ENDIF
    ENDDO

    !  Print the boundary grid data.
    write(*,*) " Boundary nodes:"
    write(*,*) "    segments = ", nbcs
    do i = 1, nbcs
       write(*,'(a9,i3,2(a11,i5))') " boundary", i, "  bnodes = ", bouns(i)%nbnode, &
            "  bfaces = ", bouns(i)%nbnode-1
    end do
    write(*,*)    

    CLOSE(1001)

  END SUBROUTINE dataInput_mpfem
!********************************************************************************

!********************************************************************************
! create outpt file in vtk legacy format
  SUBROUTINE vtkIO_legacy(fname)

    USE constants
    USE mesh_data, only : nodes,elems,wprim

    IMPLICIT NONE

    !Inputs:
    CHARACTER (LEN=80),INTENT(IN) :: fname

    !Local variables:
    INTEGER :: i,j,k
    INTEGER,PARAMETER :: etype = 5 !element type
    INTEGER :: etotal 
    CHARACTER (LEN=30) :: datasetName,fldName,meshName,endflds

    etotal = 4*nelem

    OPEN(1,FILE = fname)

    WRITE(1,"(A26)") '# vtk DataFile Version 4.2'
    WRITE(1,"(A20)") 'Incompressible Euler'
    WRITE(1,"(A5)") 'ASCII'
    WRITE(1,*) ''
    WRITE(1,"(A25)") 'DATASET UNSTRUCTURED_GRID'
    WRITE(1,"(A6,1X,I5,1X,A6)") 'POINTS',nnode,'DOUBLE'

    DO i=1,nnode
       WRITE(1,*) nodes(i)%x,nodes(i)%y,0.0d0
    ENDDO   
    WRITE(1,*) ''

    WRITE(1,"(A5,1X,I5,1X,I5)") 'CELLS',nelem,etotal
    !subtract 1 from index for vtk
    DO i=1,nelem
       WRITE(1,*) elems(i)%nnoel,elems(i)%lnode(1)-1,elems(i)%lnode(2)-1,elems(i)%lnode(3)-1 
    ENDDO   
    WRITE(1,*)''

    WRITE(1,"(A10,1X,I5)") 'CELL_TYPES',nelem
    DO i=1,nelem
       WRITE(1,"(I1)") etype
    ENDDO   
    WRITE(1,*)''

    WRITE(1,"(A10,1X,I5)") 'POINT_DATA',nnode
    WRITE(1,"(A23)") 'SCALARS pressure DOUBLE'    
    WRITE(1,"(A20)") 'LOOKUP_TABLE default'    
!!$    DO i=1,nnode
!!$       WRITE(1,*) wprim(i)%pg
!!$    ENDDO
!!$
!!$    WRITE(1,"(A23)") 'SCALARS density DOUBLE'    
    DO i=1,nnode
       WRITE(1,*) wprim(i)%d
    ENDDO
   
    WRITE(1,"(A23)") 'VECTORS velocity DOUBLE'
    DO i=1,nnode
       WRITE(1,*) wprim(i)%vx,wprim(i)%vy,0.0d0
    ENDDO   
    WRITE(1,*) ''

    CLOSE(1)
    
  END SUBROUTINE vtkIO_legacy   
!********************************************************************************

!********************************************************************************
!* Print primative variables
  subroutine printPrim(istart,iend)

    USE mesh_data, only : wprim

    IMPLICIT NONE
    
    INTEGER,INTENT(IN)                :: istart,iend !starting and ending node
    INTEGER                           :: i,j,k
    DOUBLE PRECISION                  :: KE !kinetic energy

    do i=istart,iend
       write(*,*) 'node =',i,'d =',wprim(i)%d,'vx =',wprim(i)%vx,'vy =',wprim(i)%vy,'pg =',wprim(i)%pg
    enddo

  end subroutine printPrim
!********************************************************************************

!********************************************************************************
!* Print conservative variables
  subroutine printCons(istart,iend)

    USE mesh_data, only : ucons

    IMPLICIT NONE
    
    INTEGER,INTENT(IN)                :: istart,iend !starting and ending node
    INTEGER                           :: i,j,k
    DOUBLE PRECISION                  :: KE !kinetic energy

    do i=istart,iend
       write(*,*) 'node =',i,'d =',ucons(i)%d,'mx =',ucons(i)%mx,'my =',ucons(i)%my,'en =',ucons(i)%en
    enddo

  end subroutine printCons
!********************************************************************************

END MODULE dataIO
