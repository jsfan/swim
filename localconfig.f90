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

!> localconfig.f90
!> Sharing data between main routine and signal handler
!! @authors Christian Lerrahn
module localconfig

  ! configuration options shared with swim.f90
  
  integer :: outint !< Interval for output to NetCDF file (in steps)
  integer :: outtime !< Interval for output to NetCDF file (in seconds)
  integer :: maxsteps !< Maximum number of steps to run for
  integer :: resumeint !< Interval for resume files
  integer :: finit !< flag for init/resume

  real :: maxtime !< maximum model time 
  real :: relax !< relaxation parameter for SOR solver

  character*300 :: outfile !< Name and path of output file
  character*300 :: infile !< Name and path of input (initialisation or resume) file


end module localconfig
