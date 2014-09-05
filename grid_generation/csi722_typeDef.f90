!!$ Program: csi786_typeDef.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 10/22/2011
!!$ Description: Data structure for faces


MODULE facedef
  
  INTEGER,PARAMETER :: gtype = 0  !grid type: 0:nonuniform; 1:uniform
  INTEGER,PARAMETER :: ndim = 2  !number of dimensions
  INTEGER,PARAMETER :: nfael = 3 !number of faces per element
  INTEGER,PARAMETER :: nnofa = 2 !number of nodes per face
  INTEGER,PARAMETER :: nelem_max = 1000 !maximum number of elements
  INTEGER,PARAMETER :: nface_max = 1000 !maximum number of faces
  INTEGER,PARAMETER :: npoin_max = 10000 !maximum number of points
  DOUBLE PRECISION,PARAMETER :: pi = (ATAN(1.0d0))*4.0d0 !define irrational number, pi
  DOUBLE PRECISION,PARAMETER :: tol = 1.0E-12 !define tolerance
!!$  DOUBLE PRECISION,PARAMETER :: he = 0.250d0 !height of uniform trinangular element
!!$  DOUBLE PRECISION,PARAMETER :: he2 = he*he !area of uniform trinangular element
  DOUBLE PRECISION,PARAMETER :: leps = 1.0E-6 !small value to adjust face length
  DOUBLE PRECISION,PARAMETER :: sin60 = DSIN(pi/3.0d0)
  DOUBLE PRECISION,PARAMETER :: tan60 = DTAN(pi/3.0d0)
  DOUBLE PRECISION :: he  !height of trinangular element
  DOUBLE PRECISION :: he2  !area of trinangular element


  TYPE FACE
     SEQUENCE

     INTEGER :: fid        !face id
     INTEGER :: eleid      !element id
     INTEGER :: state        !Is face part of active front? 0:no; 1:yes 
     DOUBLE PRECISION :: flen                 !artifical face length
     INTEGER,DIMENSION(nnofa) :: fnode !nodes of face
     DOUBLE PRECISION,DIMENSION(ndim) :: nhat !unit normal (inward)

  END TYPE FACE


END MODULE facedef
