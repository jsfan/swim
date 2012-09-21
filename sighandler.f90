!> SWiM - a semi-Lagrangian, semi-implicit shallow water model in
!! Cartesian coordiates
!! Copyright (C) 2008-2012 Christian Lerrahn
!!
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.

!> sighandler.f90
!> Handles system signals.
!! Currently only SIGINT (Ctrl-C) is caught.
!! @authors Christian Lerrahn
module sighandler

  use output
  use localconfig

  implicit none

contains

  !> Initialises the signal handler
  !! registering a handler for SIGINT
  subroutine sighandler_init()
    
    integer :: watchsignal
    integer :: status


    ! set up watch for following signal types
    status = watchsignal(2) ! SIGINT (Ctrl-C)

  end subroutine sighandler_init

  !> Check if any signals have been caught.
  !! @param counter Number to use for resume file written before shutdown.
  subroutine check_signal(counter)

    integer :: getlastsignal
    integer :: csig
    integer :: counter


    ! read last signal sent
    csig = getlastsignal()

    if (csig.eq.2) then ! SIGINT
       write(*,*) 'Caught interrupt signal. Saving data and exiting gracefully.'

       call write_resume_file(counter)
       call close_output_file()

       ! exit
       stop

    end if
    
  end subroutine check_signal

end module sighandler
