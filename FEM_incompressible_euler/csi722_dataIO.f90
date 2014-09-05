!!$ Program: csi722_dataIO.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 12/15/2011
!!$ Description: Module to handle input and write output in zfem format.

MODULE handleData
USE facedef

  INTERFACE dataInput
     MODULE PROCEDURE dataInput_lohner
  END INTERFACE

  INTERFACE vtkIO
     MODULE PROCEDURE vtkIO_legacy!,vtkIO_xml
  END INTERFACE

!!$  INTERFACE dataOutput
  INTERFACE zfemIO
!!$     MODULE PROCEDURE dataOutput_vecfld, dataOutput_scalfld, dataOutput_VecStruc
     MODULE PROCEDURE zfemIO_vecfld, zfemIO_scalfld, zfemIO_VecStruc
  END INTERFACE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! dataInput_lohner !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE dataInput_lohner(ndim,npoin,nnode,nelem,nboun,lnode,bnode,xyz,initVal,fname)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,npoin,nnode,nelem,nboun
    INTEGER :: i,j,k
    INTEGER :: bufferI          !buffer for integer values
    INTEGER,INTENT(OUT),DIMENSION(nnode,nelem) :: lnode
    INTEGER,INTENT(OUT),DIMENSION(2,nboun) :: bnode
    DOUBLE PRECISION,INTENT(OUT),DIMENSION(ndim,npoin) :: xyz
    DOUBLE PRECISION,INTENT(OUT),DIMENSION(4,npoin) :: initVal
    CHARACTER (LEN=50),INTENT(IN) :: fname
    CHARACTER (LEN=30) :: datasetName,fldName,meshName,endflds
    CHARACTER (LEN=50):: bufferC ! Buffer for information that is disregarded
    
    !initialize data
    lnode(:,:) = 0
    bnode(:,:) = 0
    xyz(:,:) = 0.0d0
        
    ! Open file containing input data
    OPEN(1001,FILE = fname)

    READ(1001,*) bufferI
    READ(1001,*) bufferC
    READ(1001,*) bufferC
    READ(1001,*) bufferC
    READ(1001,*) bufferC
    READ(1001,*) bufferI
    READ(1001,*) bufferC
    READ(1001,*) bufferI
    READ(1001,*) bufferC

    !read in connectivity array
    DO i=1,nelem
       READ(1001,*) bufferI,lnode(1,i),lnode(2,i),lnode(3,i)
    ENDDO

    READ(1001,*) bufferC

    !read in point coordinates
    DO i=1,npoin
       READ(1001,*) bufferI,xyz(1,i),xyz(2,i)
    ENDDO

    READ(1001,*) bufferC

    !read initial values
    DO i=1,npoin
       READ(1001,*) bufferI,initVal(1,i),initVal(2,i),initVal(3,i),initVal(4,i)
    ENDDO

    READ(1001,*) bufferC

    !read boundary points
    !bnode(1,i) is field point
    !bnode(2,i) is boundary type, 0:wall, 4:farfield
    DO i=1,nboun
       READ(1001,*) bnode(1,i),bnode(2,i)
    ENDDO

    CLOSE(1001)

  END SUBROUTINE dataInput_lohner

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! vtkIO_legacy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE vtkIO_legacy(npoin,nnode,nelem,lnode,xyz,Pr,vfld,fname)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: npoin,nnode,nelem
    INTEGER :: i,j,k
    INTEGER,PARAMETER :: etype = 5 !element type
    INTEGER :: etotal 
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    DOUBLE PRECISION,INTENT(IN),DIMENSION(3,npoin) :: xyz
    DOUBLE PRECISION,INTENT(IN),DIMENSION(npoin) :: Pr
    DOUBLE PRECISION,INTENT(IN),DIMENSION(3,npoin) :: vfld
    CHARACTER (LEN=30),INTENT(IN) :: fname
    CHARACTER (LEN=30) :: datasetName,fldName,meshName,endflds

    etotal = 4*nelem

    OPEN(1,FILE = fname)

    WRITE(1,"(A26)") '# vtk DataFile Version 4.2'
    WRITE(1,"(A20)") 'Incompressible Euler'
    WRITE(1,"(A5)") 'ASCII'
    WRITE(1,*) ''
    WRITE(1,"(A25)") 'DATASET UNSTRUCTURED_GRID'
    WRITE(1,"(A6,1X,I5,1X,A6)") 'POINTS',npoin,'DOUBLE'

    DO i=1,npoin
       WRITE(1,*) xyz(1,i),xyz(2,i),xyz(3,i)
    ENDDO   
    WRITE(1,*) ''

    WRITE(1,"(A5,1X,I5,1X,I5)") 'CELLS',nelem,etotal
    !subtract 1 from index for vtk
    DO i=1,nelem
       WRITE(1,*) nnode,lnode(1,i)-1,lnode(2,i)-1,lnode(3,i)-1 
    ENDDO   
    WRITE(1,*)''

    WRITE(1,"(A10,1X,I5)") 'CELL_TYPES',nelem
    DO i=1,nelem
       WRITE(1,"(I1)") etype
    ENDDO   
    WRITE(1,*)''

    WRITE(1,"(A10,1X,I5)") 'POINT_DATA',npoin    
    WRITE(1,"(A23)") 'SCALARS pressure DOUBLE'    
    WRITE(1,"(A20)") 'LOOKUP_TABLE default'    
    DO i=1,npoin
       WRITE(1,*) Pr(i)
    ENDDO
   
    WRITE(1,"(A23)") 'VECTORS velocity DOUBLE'
    DO i=1,npoin
       WRITE(1,*) vfld(1,i),vfld(2,i),vfld(3,i)
    ENDDO   
    WRITE(1,*) ''

    CLOSE(1)
    
  END SUBROUTINE vtkIO_legacy   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! vtkIO_xml !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  SUBROUTINE vtkIO_xml(npoin,nnode,nelem,lnode,xyz,fld,fname)
!!$    IMPLICIT NONE
!!$
!!$    INTEGER,INTENT(IN) :: npoin,nnode,nelem
!!$    INTEGER :: i,j,k
!!$    INTEGER,PARAMETER :: etype = 5 !element type
!!$    INTEGER :: etotal 
!!$    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
!!$    DOUBLE PRECISION,INTENT(IN),DIMENSION(3,npoin) :: xyz
!!$    DOUBLE PRECISION,INTENT(IN),DIMENSION(npoin) :: fld
!!$    CHARACTER (LEN=30),INTENT(IN) :: fname
!!$    CHARACTER (LEN=30) :: datasetName,fldName,meshName,endflds
!!$
!!$    etotal = 4*nelem
!!$
!!$    OPEN(1,FILE = fname)
!!$
!!$
!!$    WRITE(1,"(A26)") '# vtk DataFile Version 4.2'
!!$    WRITE(1,"(A20)") 'Incompressible Euler'
!!$    WRITE(1,"(A5)") 'ASCII'
!!$    WRITE(1,*) ''
!!$    WRITE(1,"(A25)") 'DATASET UNSTRUCTURED_GRID'
!!$    WRITE(1,"(A6,1X,I5,1X,A6)") 'POINTS',npoin,'DOUBLE'
!!$
!!$    DO i=1,npoin
!!$       WRITE(1,*) xyz(1,i),xyz(2,i),xyz(3,i)
!!$    ENDDO   
!!$    WRITE(1,*) ''
!!$
!!$    WRITE(1,"(A5,1X,I5,1X,I5)") 'CELLS',nelem,etotal
!!$    !subtract 1 from index for vtk
!!$    DO i=1,nelem
!!$       WRITE(1,*) nnode,lnode(1,i)-1,lnode(2,i)-1,lnode(3,i)-1 
!!$    ENDDO   
!!$    WRITE(1,*)''
!!$
!!$    WRITE(1,"(A10,1X,I5)") 'CELL_TYPES',nelem
!!$    DO i=1,nelem
!!$       WRITE(1,"(I1)") etype
!!$    ENDDO   
!!$    WRITE(1,*)''
!!$
!!$    WRITE(1,"(A10,1X,I5)") 'POINT_DATA',npoin    
!!$    WRITE(1,"(A23)") 'SCALARS pressure DOUBLE'    
!!$    WRITE(1,"(A20)") 'LOOKUP_TABLE default'    
!!$    DO i=1,npoin
!!$       WRITE(1,*) fld(i)
!!$    ENDDO   
!!$    WRITE(1,*) ''
!!$
!!$    CLOSE(1)
!!$
!!$  END SUBROUTINE vtkIO_xml   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! zfemIO_scalfld !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE zfemIO_scalfld(ndim,npoin,nnode,nelem,lnode,xyz,fld,fname)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,npoin,nnode,nelem
    INTEGER :: i,j,k
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    DOUBLE PRECISION,INTENT(IN),DIMENSION(3,npoin) :: xyz
    DOUBLE PRECISION,INTENT(IN),DIMENSION(npoin) :: fld
    CHARACTER (LEN=30),INTENT(IN) :: fname
    CHARACTER (LEN=30) :: datasetName,fldName,meshName,endflds

    WRITE(datasetName,11)
11 FORMAT('ptstri')
    WRITE(fldName,12)
12 FORMAT('fldtri')
    WRITE(meshName,13)
13 FORMAT('elemtri')
    WRITE(endflds,14)
14 FORMAT('END_FIELDS')

    OPEN(1,FILE = fname)


    WRITE(1,*)'CONTINUOUS'
    WRITE(1,*),datasetName
    WRITE(1,*),'POINTS'
    WRITE(1,*),npoin

    DO i=1,npoin
       WRITE(1,*),xyz(1,i),xyz(2,i),xyz(3,i)
    ENDDO   

    WRITE(1,*),'REAL'
    WRITE(1,*),fldName
    WRITE(1,*),ndim

    DO i=1,npoin
       WRITE(1,*),fld(i)
    ENDDO

    WRITE(1,*),endflds
    WRITE(1,*),'TRIANGLE'
    WRITE(1,*),meshName
    WRITE(1,*),nelem

    DO i=1,nelem
       WRITE(1,*),lnode(1,i),lnode(2,i),lnode(3,i) 
    ENDDO   
    WRITE(1,*),'END'

    CLOSE(1)

  END SUBROUTINE zfemIO_scalfld   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! zfemIO_vecfld !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE zfemIO_vecfld(ndim,npoin,nnode,nelem,lnode,xyz,vfld,fname)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,npoin,nnode,nelem
    INTEGER :: i,j,k
    INTEGER,INTENT(IN),DIMENSION(nnode,nelem) :: lnode
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: xyz
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,npoin) :: vfld
    CHARACTER (LEN=30),INTENT(IN) :: fname
    CHARACTER (LEN=30) :: datasetName,fldName,meshName,endflds

    WRITE(datasetName,11)
11 FORMAT('ptstri')
    WRITE(fldName,12)
12 FORMAT('gradtri')
    WRITE(meshName,13)
13 FORMAT('elemtri')
    WRITE(endflds,14)
14 FORMAT('END_FIELDS')

    OPEN(1,FILE = fname)

    WRITE(1,*),'CONTINUOUS'
    WRITE(1,*),datasetName
    WRITE(1,*),'POINTS'
    WRITE(1,*),npoin

    DO i=1,npoin
       WRITE(1,*),xyz(1,i),xyz(2,i),xyz(3,i)
    ENDDO   

    WRITE(1,*),'REAL'
    WRITE(1,*),fldName
    WRITE(1,*),ndim

    DO i=1,npoin
       WRITE(1,*),vfld(1,i),vfld(2,i),vfld(3,i) 
    ENDDO

    WRITE(1,*),endflds
    WRITE(1,*),'TRIANGLE'
    WRITE(1,*),meshName
    WRITE(1,*),nelem

    DO i=1,nelem
       WRITE(1,*),lnode(1,i),lnode(2,i),lnode(3,i) 
    ENDDO   
    WRITE(1,*),'END'

    CLOSE(1)

  END SUBROUTINE zfemIO_vecfld   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! zfemIO_VecStruc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE zfemIO_VecStruc(ndim,nx,ny,nz,nt,xyz,vfld,fname)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim,nx,ny,nz,nt
    INTEGER :: i,j,k
    INTEGER:: ndxyz,npx,npy,npz
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,nx) :: xyz
    DOUBLE PRECISION,INTENT(IN),DIMENSION(ndim,nt) :: vfld
    CHARACTER (LEN=30),INTENT(IN) :: fname
    CHARACTER (LEN=30) :: datasetName,fldName,meshName,endflds

    ndxyz = 0
    npx = nx
    npy = ny
    npz = nz
    IF(nx .NE. 1)THEN
       ndxyz = ndxyz + 1
    ELSE
       npx = 0
    ENDIF   

    IF(ny .NE. 1)THEN
       ndxyz = ndxyz + 1
    ELSE
       npy = 0
    ENDIF   

    IF(nz .NE. 1)THEN
       ndxyz = ndxyz + 1
    ELSE
       npz = 0
    ENDIF   

    WRITE(datasetName,11)
11 FORMAT('ptsstruc')
    WRITE(fldName,12)
12 FORMAT('gradstruc')
    WRITE(meshName,13)
13 FORMAT('elemstruc')
    WRITE(endflds,14)
14 FORMAT('END_FIELDS')

    OPEN(1,FILE = fname)

    WRITE(1,*),'CONTINUOUS'
    WRITE(1,*),datasetName
    WRITE(1,*),'POINTS'
    WRITE(1,*),nt

    DO i=1,nx
       DO j=1,ny
          WRITE(1,*),xyz(1,i),xyz(2,j),xyz(3,1)
       ENDDO
    ENDDO   

    WRITE(1,*),'REAL'
    WRITE(1,*),fldName
    WRITE(1,*),ndim

    DO i=1,nt
       WRITE(1,*),vfld(1,i),vfld(2,i),vfld(3,i) 
    ENDDO

    WRITE(1,*),endflds
    WRITE(1,*),'STRUCTURED'
    WRITE(1,*),meshName
    WRITE(1,*),ndxyz,npx,npy,npz

    WRITE(1,*),'END'

    CLOSE(1)

  END SUBROUTINE zfemIO_VecStruc

END MODULE handleData
