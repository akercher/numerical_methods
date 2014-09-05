!!$ Program: csi722_typeDef.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 12/15/2011
!!$ Description: data structure for faces


MODULE facedef
  
  INTEGER,PARAMETER :: nnofa = 2 !number of nodes per face
  INTEGER,PARAMETER :: nfael = 2 !number of faces per element
  DOUBLE PRECISION,PARAMETER :: pi = (ATAN(1.0d0))*4.0d0
  DOUBLE PRECISION,PARAMETER :: gamma = 1.40d0 !ratio of specific heat, Cp/Cv 
  DOUBLE PRECISION,PARAMETER :: Minf = 1.0d0 !far-field Mach Number
  DOUBLE PRECISION,PARAMETER :: alpha = 5.0d0*pi/180.0d0 !angle of attack in radians
  DOUBLE PRECISION,PARAMETER :: Prinf = 1.0d0 
  DOUBLE PRECISION,PARAMETER :: cinf = 1.0d0 
  DOUBLE PRECISION,PARAMETER :: rhoinf = gamma*Prinf/(cinf*cinf)
  DOUBLE PRECISION,PARAMETER :: vinf = Minf*cinf 
  DOUBLE PRECISION,PARAMETER :: vxinf = vinf*DCOS(alpha)
  DOUBLE PRECISION,PARAMETER :: vyinf = vinf*DSIN(alpha)

  TYPE FACE
     SEQUENCE

     INTEGER :: eleid           !element id
     INTEGER :: btype           !type of boundary point (0:wall,4:farfield)
     INTEGER :: flow            !flow type (0:inflow,1:outflow,-1:wall,-2:wingtip of wall)
     DOUBLE PRECISION :: vninf  !normal value at infinity
     DOUBLE PRECISION :: vtinf  !tangent value at infinity
     INTEGER,DIMENSION(nnofa) :: fnode !nodes of face
     DOUBLE PRECISION,DIMENSION(2) :: nhat !unit normal (outward)
     DOUBLE PRECISION,DIMENSION(2) :: that !unit tangent (conterclockwise)

  END TYPE FACE

  TYPE NODE
     SEQUENCE

     INTEGER :: bpt             !boundary point
     INTEGER :: btype           !type of boundary point (0:wall,4:farfield)
     INTEGER :: flow            !flow type (0:inflow,1:outflow)
     DOUBLE PRECISION :: vninf  !normal value at infinity
     DOUBLE PRECISION :: vtinf  !tangent value at infinity
     INTEGER,DIMENSION(2) :: fid           !face ids
     DOUBLE PRECISION,DIMENSION(2) :: nhat !unit normal (outward)
     DOUBLE PRECISION,DIMENSION(2) :: that !unit tangent (conterclockwise)

  END TYPE NODE

CONTAINS


END MODULE facedef
