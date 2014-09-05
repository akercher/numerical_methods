!!$ Program: Euler2D_timeCLF.f90
!!$ Author: Andrew Kercher
!!$ Last Updated: 02/09/2012
!!$ Description: Module includes subroutines for local time step calculation.
!!$ Acknowledgment: The following algorithms and diagrams were taken from 
!!$                 or heavily based off of code originally written by 
!!$                 Dr. Katate Masatsuka (info[at]cfdbooks.com).


MODULE timeCLF

private                        

public :: local_time_step

CONTAINS

  subroutine local_time_step(dtg)

    use constants
    use mesh_data, only : nodes

    IMPLICIT NONE

    !Local variables
    integer  :: i,j,k
    DOUBLE PRECISION,INTENT(OUT) :: dtg
    
    dtg = 1.0E+16
    do i = 1,nnode

       ! Local time step: dt = volume/sum(0.5*max_wave_speed*face_area).
       nodes(i)%dt = nodes(i)%ml / nodes(i)%nwsd

       ! Keep the minimum dtl
       dtg = min(dtg,nodes(i)%dt)

    end do

  end subroutine local_time_step

END MODULE timeCLF
